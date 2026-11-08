"""Count target-word sets and exact targets in paired reads.

Each detector is built for one target. Evaluation state belongs to one call,
so repeated simulations cannot leak evidence. A word set does not encode word
order or multiplicity; exact full-target matching explicitly requires both.
"""
from collections import deque

RC = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')
DNA = frozenset('ACGTRYSWKMBDHVN')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def rc(sequence):
    return sequence.translate(RC)[::-1]


def canonical(sequence):
    return min(sequence, rc(sequence))


class LiteralMatcher:
    """Aho-Corasick literal matches, including overlaps and both orientations.

    Canonical patterns have one identifier even when palindromic. Ambiguous
    read bases reset matching. The Q20 flag requires every base in the whole
    matched interval to have Phred+33 >= 20.
    """
    def __init__(self, patterns):
        self.patterns = tuple(sorted(set(patterns)))
        require(self.patterns and all(p and set(p) <= set('ACGT') and p == canonical(p)
                                     for p in self.patterns),
                'Patterns must be nonempty canonical ACGT strings')
        edges, fail, outputs = [{}], [0], [[]]
        for pid, pattern in enumerate(self.patterns):
            for word in sorted({pattern, rc(pattern)}):
                node = 0
                for base in word:
                    if base not in edges[node]:
                        edges[node][base] = len(edges)
                        edges.append({})
                        fail.append(0)
                        outputs.append([])
                    node = edges[node][base]
                outputs[node].append(pid)
        pending = deque(edges[0].values())
        while pending:
            node = pending.popleft()
            for base, child in edges[node].items():
                pending.append(child)
                fallback = fail[node]
                while fallback and base not in edges[fallback]:
                    fallback = fail[fallback]
                fail[child] = edges[fallback].get(base, 0)
                outputs[child].extend(outputs[fail[child]])
        self.edges = tuple(edges)
        self.fail = tuple(fail)
        self.outputs = tuple(tuple(o) for o in outputs)
        self.ids = {p: i for i, p in enumerate(self.patterns)}
        self.lengths = tuple(map(len, self.patterns))

    def scan(self, sequence, quality):
        require(isinstance(sequence, str) and sequence and set(sequence) <= DNA,
                'Read sequence must be a nonempty uppercase IUPAC string')
        require(isinstance(quality, str) and len(sequence) == len(quality)
                and all(33 <= ord(q) <= 126 for q in quality),
                'Quality must be equal-length printable Phred+33')
        low = [0]
        for q in quality:
            low.append(low[-1] + (ord(q) < 53))
        node = 0
        hits = []
        for end, base in enumerate(sequence, 1):
            while node and base not in self.edges[node]:
                node = self.fail[node]
            node = self.edges[node].get(base, 0)
            for pid in self.outputs[node]:
                start = end - self.lengths[pid]
                hits.append((pid, start, end, low[end] == low[start]))
        return hits


class Detector:
    """Evaluate the identical paired-read observations with all criteria.

    Pair counts count observations, at most once per word or criterion per
    pair. Signature counts collapse identical ordered (R1 sequence, R2
    sequence), ignoring qualities. A same-pair/read/exact criterion must pass
    in an individual observation before its signature is credited; quality
    evidence from duplicate observations is never spliced for those criteria.
    Pooled word signatures, by definition, can accrue across observations.
    """
    def __init__(self, sequence, k_values=(31, 41, 51, 61, 71, 81)):
        require(isinstance(sequence, str) and sequence and set(sequence) <= set('ACGT'),
                'Target must be nonempty uppercase ACGT')
        values = tuple(k_values)
        require(all(isinstance(k, int) and not isinstance(k, bool) and k > 0 for k in values),
                'k values must be positive integers')
        self.sequence = sequence
        self.k_values = tuple(sorted({k for k in values if k <= len(sequence)} | {len(sequence)}))
        self.words = tuple(tuple(sorted({canonical(sequence[i:i+k])
                                        for i in range(len(sequence)-k+1)}))
                           for k in self.k_values)
        self.matcher = LiteralMatcher(word for words in self.words for word in words)
        self.pattern_ids = tuple(tuple(self.matcher.ids[word] for word in words) for words in self.words)
        owners = [None] * len(self.matcher.patterns)
        for ki, ids in enumerate(self.pattern_ids):
            for bit, pid in enumerate(ids):
                require(owners[pid] is None, 'Internal word ownership collision')
                owners[pid] = (ki, 1 << bit)
        self.owners = tuple(owners)
        self.full_masks = tuple((1 << len(words)) - 1 for words in self.words)
        self.exact_pid = self.matcher.ids[canonical(sequence)]

    def evaluate(self, reads):
        n_patterns, n_k = len(self.matcher.patterns), len(self.k_values)
        pair_counts = [[0] * n_patterns for _ in range(2)]
        word_signatures = [[set() for _ in range(n_patterns)] for _ in range(2)]
        same_pair_counts = [[0] * n_k for _ in range(2)]
        same_read_counts = [[0] * n_k for _ in range(2)]
        same_pair_signatures = [[set() for _ in range(n_k)] for _ in range(2)]
        same_read_signatures = [[set() for _ in range(n_k)] for _ in range(2)]
        signatures = {}
        for row in reads:
            seqs = (row['r1_sequence'], row['r2_sequence'])
            quals = (row['r1_quality'], row['r2_quality'])
            hits = [self.matcher.scan(s, q) for s, q in zip(seqs, quals)]
            sid = signatures.setdefault(seqs, len(signatures))
            for qi in range(2):
                seen = [{pid for pid, start, end, high in hh if qi == 0 or high} for hh in hits]
                masks = [[0] * n_k for _ in range(2)]
                for mate in range(2):
                    for pid in seen[mate]:
                        ki, bit = self.owners[pid]
                        masks[mate][ki] |= bit
                for pid in seen[0] | seen[1]:
                    pair_counts[qi][pid] += 1
                    word_signatures[qi][pid].add(sid)
                for ki, full in enumerate(self.full_masks):
                    left, right = masks[0][ki], masks[1][ki]
                    if (left | right) == full:
                        same_pair_counts[qi][ki] += 1
                        same_pair_signatures[qi][ki].add(sid)
                    if left == full or right == full:
                        same_read_counts[qi][ki] += 1
                        same_read_signatures[qi][ki].add(sid)
        result = []
        for ki, k in enumerate(self.k_values):
            ids = self.pattern_ids[ki]
            for qi, quality in enumerate(('raw', 'q20')):
                counts = dict(
                    minimum_word_pairs=min(pair_counts[qi][pid] for pid in ids),
                    minimum_word_signatures=min(len(word_signatures[qi][pid]) for pid in ids),
                    same_pair_pairs=same_pair_counts[qi][ki],
                    same_pair_signatures=len(same_pair_signatures[qi][ki]),
                    same_read_pairs=same_read_counts[qi][ki],
                    same_read_signatures=len(same_read_signatures[qi][ki]),
                    exact_pairs=pair_counts[qi][self.exact_pid],
                    exact_signatures=len(word_signatures[qi][self.exact_pid]))
                for unit in ('pairs', 'signatures'):
                    require(counts['exact_' + unit] <= counts['same_read_' + unit]
                            <= counts['same_pair_' + unit] <= counts['minimum_word_' + unit],
                            'Evidence nesting invariant failed')
                    if k == len(self.sequence):
                        require(len({counts[name + '_' + unit] for name in
                                     ('minimum_word', 'same_pair', 'same_read', 'exact')}) == 1,
                                'Full-length word and exact match must be equivalent')
                result.append(dict(k=k, span_length=len(self.sequence), quality=quality,
                                   distinct_canonical_words=len(ids), **counts))
        return result
