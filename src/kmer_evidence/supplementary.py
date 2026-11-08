"""Supplementary catalog analyses computed from the supplied derived counters.

Fraction rule
    Sequence indexes often call a query present when a chosen fraction of its
    k-mers is found. Here a target passes the fraction rule when at least the
    stated fraction of its distinct canonical k-mers each reach the support
    threshold. A required fraction of 1.0 is the pooled rule of the main
    comparison (the minimum k-mer support reaches the threshold).

Target repeat profile
    The longest substring that occurs more than once in a target, counting
    both strands. Same-read and exact counts can differ only when a read holds
    every target k-mer without holding the contiguous target; internal repeats
    are the simplest way for that to happen.

Selected-case structure
    How the selected pooled-only cases are distributed over target states and
    anchor pairs, using the supplied case ledger.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
import csv
import gzip
import json
from pathlib import Path

FRACTIONS = (0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
_COMP = str.maketrans('ACGT', 'TGCA')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def _rows(path):
    with gzip.open(path, 'rt', newline='') as handle:
        yield from csv.DictReader(handle)


def reverse_complement(sequence):
    return sequence.translate(_COMP)[::-1]


def canonical(sequence):
    return min(sequence, reverse_complement(sequence))


def longest_repeat(sequence):
    """Length of the longest substring occurring at least twice on either strand."""
    require(sequence and set(sequence) <= set('ACGT'), 'Target must be an ACGT string')
    best = 0
    for length in range(1, len(sequence)):
        seen = set()
        repeated = False
        for start in range(len(sequence) - length + 1):
            word = canonical(sequence[start:start + length])
            if word in seen:
                repeated = True
                break
            seen.add(word)
        if not repeated:
            break
        best = length
    return best


@lru_cache(maxsize=32)
def _fraction_ratio(text):
    """Represent a finite fraction in [0, 1] without rounding its precision."""
    ratio = Fraction(text)
    require(0 <= ratio <= 1, 'Required fraction must be between 0 and 1')
    return ratio.numerator, ratio.denominator


def passes_fraction(supports, threshold, fraction):
    """Compare the observed share with the decimal value of `fraction` exactly."""
    require(supports, 'A target needs at least one k-mer')
    numerator, denominator = _fraction_ratio(str(fraction))
    passing = sum(value >= threshold for value in supports)
    return passing * denominator >= numerator * len(supports)


def fraction_rule_table(source, protocol, fractions=FRACTIONS):
    source = Path(source)
    k_values = [int(k) for k in protocol['k_values']]
    thresholds = [int(t) for t in protocol['support_thresholds']]
    states = {r['state_id']: r['sequence'] for r in _rows(source / 'target_sequences.csv.gz')}
    queries = {}
    for r in _rows(source / 'queries.csv.gz'):
        queries[r['pattern_id']] = r['sequence']
    groups = {}
    labels = {}
    for r in _rows(source / 'groups.csv.gz'):
        words = json.loads(r['word_ids'])
        seq = states[r['state_id']]
        k = int(r['k'])
        require({queries[w] for w in words} == {canonical(seq[i:i + k]) for i in range(len(seq) - k + 1)},
                'Group k-mer set differs from its target: ' + r['group_id'])
        groups[r['group_id']] = (r['state_id'], words, r['span_query_id'])
        labels[r['group_id']] = json.loads(r['baseline_labels'])
    support = {}
    for r in _rows(source / 'word_counts.csv.gz'):
        support[r['pattern_id'], r['pi_accession']] = (int(r['raw_pair_support']), int(r['q20_pair_support']))
    accessions = sorted({pi for _, pi in support})
    require(len(accessions) == protocol['expected_accessions'], 'Accession count differs')
    cells = defaultdict(dict)
    for gid, (sid, words, span) in groups.items():
        for pi in accessions:
            for qi, quality in enumerate(('raw', 'q20')):
                values = [support[w, pi][qi] for w in words]
                exact = support[span, pi][qi]
                for label in labels[gid]:
                    cells[label][pi, sid, quality] = (values, exact)
    fixed = ['k' + str(k) for k in k_values]
    common = set.intersection(*({sid for _, sid, _ in cells[b]} for b in fixed))
    rows = []
    for scope in ('ALL_CATALOG', 'COMMON_K31_K81'):
        for baseline in fixed + ['near_full']:
            for quality in ('raw', 'q20'):
                selected = [(key, value) for key, value in cells[baseline].items()
                            if key[2] == quality and (scope == 'ALL_CATALOG' or key[1] in common)]
                for threshold in thresholds:
                    for fraction in fractions:
                        passing = without_exact = exact_positive = 0
                        for _, (values, exact) in selected:
                            ok = passes_fraction(values, threshold, fraction)
                            has_exact = exact >= threshold
                            passing += ok
                            exact_positive += has_exact
                            without_exact += ok and not has_exact
                            require(not has_exact or ok, 'Exact support without k-mer support')
                        rows.append(dict(scope=scope, baseline=baseline, quality=quality,
                                         minimum_support=threshold, required_kmer_fraction=fraction,
                                         eligible_cells=len(selected), fraction_rule_positive_cells=passing,
                                         exact_positive_cells=exact_positive,
                                         fraction_rule_positive_without_exact_cells=without_exact,
                                         count_basis='READ_PAIR_OBSERVATIONS'))
    return rows


def target_repeat_profile(source):
    rows = []
    for r in _rows(Path(source) / 'target_sequences.csv.gz'):
        seq = r['sequence']
        rows.append(dict(state_id=r['state_id'], pair_id=r['pair_id'], span_length=len(seq),
                         longest_repeated_substring_both_strands=longest_repeat(seq)))
    rows.sort(key=lambda row: row['state_id'])
    return rows


def repeat_summary(profile):
    lengths = [row['longest_repeated_substring_both_strands'] for row in profile]
    return dict(states=len(lengths), maximum=max(lengths),
                states_at_least_20=sum(x >= 20 for x in lengths),
                states_at_least_30=sum(x >= 30 for x in lengths),
                states_below_20=sum(x < 20 for x in lengths))


def selected_case_structure(ledger_path):
    ledger = list(_rows(ledger_path))
    by_label = defaultdict(list)
    for row in ledger:
        by_label[row['primary_evidence_label']].append(row)
    pair_counts = Counter(row['pair_id'] for row in ledger)
    return dict(
        cases=len(ledger), states=len({r['state_id'] for r in ledger}),
        anchor_pairs=len(pair_counts),
        largest_anchor_pair_contributions=pair_counts.most_common(3),
        labels={label: dict(cases=len(rows), states=len({r['state_id'] for r in rows}),
                            anchor_pairs=len({r['pair_id'] for r in rows}),
                            accessions=len({r['pi_accession'] for r in rows}),
                            raw_exact_pairs_zero=sum(int(r['raw_exact_pairs']) == 0 for r in rows),
                            one_q20_exact_signature=sum(r['one_q20_signature_below_two'] == 'True' for r in rows),
                            reference_placement=dict(sorted(Counter(r['placement_status'] for r in rows).items())))
                for label, rows in sorted(by_label.items())})
