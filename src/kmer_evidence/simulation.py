"""Generate fixed synthetic templates and paired reads with known template truth.

Truth is defined on the complete generating template mixture, BEFORE sampling
or sequencing error. These local 600 bp constructions are controlled stress
tests; they do not estimate genomic prevalence or biological allele frequency.
"""
import hashlib
import json
import math
import random


BASES = 'ACGT'
RC_TABLE = str.maketrans('ACGT', 'TGCA')
MODELS = ('present_target', 'present_diluted10', 'absent_background',
          'absent_single_snv', 'absent_complement_internal',
          'absent_complement_endpoints')
NOISE_MODELS = ('clean', 'error_1pct', 'quality_patch')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def rc(sequence):
    return sequence.translate(RC_TABLE)[::-1]


def canonical(sequence):
    return min(sequence, rc(sequence))


def stable_seed(*parts):
    """Type-aware deterministic seed; independent of process hash/random state."""
    payload = json.dumps(parts, ensure_ascii=True, sort_keys=True,
                         separators=(',', ':'), allow_nan=False).encode('ascii')
    return int.from_bytes(hashlib.sha256(payload).digest(), 'big')


def words(sequence, k):
    return {canonical(sequence[i:i+k]) for i in range(len(sequence)-k+1)}


def _contains(template, target):
    return target in template or rc(target) in template


def _replace(sequence, position, replacement):
    return sequence[:position] + replacement + sequence[position+1:]


def _different(base, rng):
    return rng.choice([b for b in BASES if b != base])


def _sha(value):
    return hashlib.sha256(value.encode('ascii')).hexdigest()


def validate_template_record(record):
    """Recompute truth on whole templates, including junctions and both strands."""
    target = record['target_sequence']
    templates = record['templates']
    weights = record['weights']
    require(target and set(target) <= set(BASES), 'Target must be ACGT')
    require(templates and len(templates) == len(weights), 'Invalid template weights')
    require(all(math.isfinite(w) and w > 0 for w in weights)
            and abs(sum(weights)-1.0) < 1e-12, 'Weights must be positive and sum to one')
    require(all(set(t) <= set(BASES) and len(t) == record['template_length']
                for t in templates), 'Invalid synthetic template')
    observed = [_contains(t, target) for t in templates]
    require(observed == record['template_contains_target'],
            'Template-level truth metadata differs from literal sequence')
    require(any(observed) == record['truth_present'], 'Mixture truth does not match templates')
    return observed


def build_templates(sequence, target_id, master_seed):
    """Construct all six specified models, or fail without silently dropping any.

    Internal complement variants mutate target positions floor(L/4) and
    floor(3L/4), respectively. Endpoint complement variants mutate positions
    0 and L-1. A different template supplies every unmutated k-mer; no template
    in either absent mixture contains the complete target on either strand.
    """
    require(isinstance(sequence, str) and set(sequence) <= set(BASES),
            'Target must be an unambiguous ACGT string')
    length = len(sequence)
    require(62 <= length <= 150, 'Supported target length is 62 through 150')
    target_start = (600-length)//2
    positions = (length//4, (3*length)//4)
    require(positions[1]-positions[0] >= 31, 'Internal variants cannot complement k31')
    target_words31 = words(sequence, 31)
    construction_seed = stable_seed('CODE25.templates.v1', master_seed, target_id, sequence)
    rng = random.Random(construction_seed)
    for attempt in range(1, 101):
        left = ''.join(rng.choice(BASES) for _ in range(target_start))
        right = ''.join(rng.choice(BASES) for _ in range(600-target_start-length))
        background_core = ''.join(rng.choice(BASES) for _ in range(length))
        single_site = length//2
        single = _replace(sequence, single_site, _different(sequence[single_site], rng))
        internal = [_replace(sequence, p, _different(sequence[p], rng)) for p in positions]
        endpoints = [_replace(sequence, p, _different(sequence[p], rng)) for p in (0, length-1)]
        all_cores = [sequence, background_core, single] + internal + endpoints
        all_templates = [left+core+right for core in all_cores]
        # Repetitive targets can remain present after an edit, or reappear at a
        # flank junction. Validate before reads exist and regenerate boundedly.
        if not _contains(all_templates[0], sequence):
            continue
        if any(_contains(t, sequence) for t in all_templates[1:]):
            continue
        if not target_words31 <= words(internal[0], 31) | words(internal[1], 31):
            continue
        if not all(words(sequence, k) <= words(endpoints[0], k) | words(endpoints[1], k)
                   for k in range(1, length)):
            continue
        break
    else:
        raise ValueError('Unable to construct all truth-validated controls in 100 attempts: '+str(target_id))

    specifications = {
        'present_target': ([0], [1.0], [], 'exact_target'),
        'present_diluted10': ([0, 1], [0.1, 0.9], [], 'ten_percent_target_template_weight'),
        'absent_background': ([1], [1.0], [], 'random_background'),
        'absent_single_snv': ([2], [1.0], [single_site], 'one_central_substitution'),
        'absent_complement_internal': ([3, 4], [0.5, 0.5], list(positions), 'designed_k31_complementation_stress'),
        'absent_complement_endpoints': ([5, 6], [0.5, 0.5], [0, length-1], 'designed_all_subspan_k_complementation_stress'),
    }
    records = {}
    for model in MODELS:
        indexes, weights, sites, label = specifications[model]
        templates = [all_templates[i] for i in indexes]
        contains = [_contains(t, sequence) for t in templates]
        cores = [all_cores[i] for i in indexes]
        record = dict(
            template_model=model, target_id=target_id, target_sequence=sequence,
            template_length=600, target_start0=target_start,
            target_end0=target_start+length, sampling_center0=300,
            templates=templates, weights=weights, truth_present=any(contains),
            template_contains_target=contains,
            construction=dict(
                label=label, construction_seed_hex=hex(construction_seed), attempt=attempt,
                target_relative_edit_positions0=sites, core_sequences=cores,
                core_sha256=[_sha(x) for x in cores],
                template_sha256=[_sha(x) for x in templates],
                shared_left_flank_sha256=_sha(left), shared_right_flank_sha256=_sha(right),
                target_k31_word_set_covered_by_core_union=target_words31 <= set().union(*(words(c, 31) for c in cores)),
                full_template_both_strand_truth_checked=True,
                truth_scope='complete_generating_template_mixture_before_read_sampling_or_errors'))
        validate_template_record(record)
        records[model] = record
    return records


def _choose_index(weights, rng):
    value = rng.random()
    total = 0.0
    for i, weight in enumerate(weights):
        total += weight
        if value < total:
            return i
    return len(weights)-1  # Floating point boundary, not an extra stochastic draw.


def _noise(sequence, model, rng):
    if model == 'clean':
        return sequence, 'I'*len(sequence), 0, 0
    if model == 'error_1pct':
        bases = []
        errors = 0
        for base in sequence:
            if rng.random() < 0.01:
                base = _different(base, rng)
                errors += 1
            bases.append(base)
        # Q20 corresponds to the specified 1% substitution probability. Errors
        # may accidentally restore a negative template's target: retain them.
        return ''.join(bases), '5'*len(sequence), errors, 0
    if model == 'quality_patch':
        width = min(10, len(sequence))
        start = rng.randrange(len(sequence)-width+1)
        quality = 'I'*start + '+'*width + 'I'*(len(sequence)-start-width)
        # Sequence is intentionally unmodified: this diagnostic isolates the
        # quality gate and is NOT a calibrated sequencing-error distribution.
        return sequence, quality, 0, width
    raise ValueError('Unknown noise model: '+str(model))


def simulate_cell(templates_record, seed, pair_count, read_length, noise_model):
    """Return (paired_read_dicts, audit), without conditioning on observed hits.

    A random local fragment contains the target-region center. FR mates are
    the two ends of that fragment, on a randomly chosen strand. Overlapping
    mates are permitted. This is locus-focused sampling, not whole-genome
    depth or an empirically fitted insert-size distribution.
    """
    require(isinstance(pair_count, int) and pair_count >= 1, 'pair_count must be positive')
    require(isinstance(read_length, int) and 1 <= read_length <= 280, 'Invalid read length')
    require(noise_model in NOISE_MODELS, 'Invalid noise model')
    validate_template_record(templates_record)
    rng = random.Random(seed)
    # Independent error stream avoids letting the number of sequencing errors
    # alter later fragment/template draws. The caller can also reuse cell seeds
    # across noise models for strictly matched origins if that is desired.
    noise_rng = random.Random(stable_seed('CODE25.noise.v1', seed, noise_model))
    templates = templates_record['templates']
    weights = templates_record['weights']
    center = templates_record['sampling_center0']
    min_fragment = max(180, read_length+20)
    reads = []
    read_digest = hashlib.sha256()
    origin_digest = hashlib.sha256()
    errors = low_quality = target_pairs = 0
    per_template = [0]*len(templates)
    fragments = []
    for pair_index in range(pair_count):
        template_index = _choose_index(weights, rng)
        template = templates[template_index]
        fragment_length = rng.randint(min_fragment, 300)
        lo = max(0, center-fragment_length+1)
        hi = min(center, len(template)-fragment_length)
        require(lo <= hi, 'No valid center-covering fragment start')
        start = rng.randint(lo, hi)
        end = start+fragment_length
        strand = '+' if rng.randrange(2) == 0 else '-'
        fragment = template[start:end]
        if strand == '-':
            fragment = rc(fragment)
        first, second = fragment[:read_length], rc(fragment[-read_length:])
        first, q1, e1, l1 = _noise(first, noise_model, noise_rng)
        second, q2, e2, l2 = _noise(second, noise_model, noise_rng)
        row = dict(r1_sequence=first, r2_sequence=second,
                   r1_quality=q1, r2_quality=q2,
                   origin_template_index=template_index,
                   origin_fragment_start0=start, origin_fragment_end0=end,
                   origin_strand=strand)
        reads.append(row)
        observed = [first, second, q1, q2]
        origin = [template_index, start, end, strand]
        read_digest.update((json.dumps(observed, separators=(',', ':'))+'\n').encode('ascii'))
        origin_digest.update((json.dumps(origin, separators=(',', ':'))+'\n').encode('ascii'))
        errors += e1+e2
        low_quality += l1+l2
        target_pairs += int(templates_record['template_contains_target'][template_index])
        per_template[template_index] += 1
        fragments.append(fragment_length)
    audit = dict(
        truth_present=templates_record['truth_present'], pair_count=pair_count,
        read_length=read_length, noise_model=noise_model,
        selected_target_template_pairs=target_pairs,
        selected_pairs_by_template=per_template, substitution_error_bases=errors,
        low_quality_bases=low_quality, sequenced_bases=2*pair_count*read_length,
        minimum_fragment_length=min(fragments), maximum_fragment_length=max(fragments),
        reads_and_quality_sha256=read_digest.hexdigest(), origin_sha256=origin_digest.hexdigest(),
        simulation_seed_hex=hex(seed),
        read_sampling_scope='local_fragments_covering_target_region_center_not_genome_depth',
        post_observation_truth_filtering=False)
    return reads, audit
