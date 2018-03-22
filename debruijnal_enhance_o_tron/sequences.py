import random

complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 


def comp(base):
    base.upper()
    return complements.get(base, base)


def revcomp(sequence):
    sequence.upper()
    bases = reversed([complements.get(base, base) for base in sequence])
    return ''.join(bases)


def mutate_base(base):
    if base in 'AT':
        return random.choice('GC')
    elif base in 'GC':
        return random.choice('AT')
    else:
        assert False, 'bad base'


def mutate_sequence(sequence, N=1):
    sequence = list(sequence)
    positions = random.sample(range(len(sequence)), N)

    for i in positions:
        sequence[i] = mutate_base(sequence[i])

    return ''.join(sequence)


def mutate_position(sequence, pos):
    sequence = list(sequence)
    sequence[pos] = mutate_base(sequence[pos])
    return ''.join(sequence)


def get_random_sequence(length, ksize, exclude=None, seen=None):
    '''Generate a random (non-looping) nucleotide sequence.

    To be non-overlapping, the sequence should not include any repeated
    length K-1 k-mers.

    Args:
        exclude (str): If not None, add the k-mers from this sequence to the
        seen set.
        seen (set): An existing set of seen k-mers.

    Returns:
        str: A random non-looping sequence.
    '''

    seen = set() if seen is None else seen.copy()

    def add_seen(kmer):
        seen.add(kmer)
        seen.add(revcomp(kmer))

    if exclude is not None:
        for pos in range(0, len(exclude) - ksize):
            add_seen(exclude[pos:pos + ksize - 1])


    seq = [random.choice('ACGT') for _ in range(ksize - 1)]  # do first K-1 bases
    add_seen(''.join(seq))

    i = 0
    while(len(seq) < length):
        if i > length:
            raise ValueError('K too small for request length.')
        i += 1
        next_base = random.choice('ACGT')
        next_kmer = ''.join(seq[-ksize + 2:] + [next_base])
        assert len(next_kmer) == ksize - 1
        if (next_kmer) not in seen:
            seq.append(next_base)
            add_seen(next_kmer)
        else:
            continue

    return ''.join(seq)


def reads(sequence, ksize, L=100, N=100, dbg_cover=False):
    '''Yields simulated reads over the given sequence.
    '''
    positions = list(range(len(sequence) - L))
    if dbg_cover is True:
        for start in range(0, len(sequence), ksize):
            read = sequence[start:start + L]
            if len(read) < ksize:
                read = sequence[-L:]
            yield read
            N -= 1
    if N < 0:
        return
    for i in range(N):
        start = random.choice(positions)
        yield sequence[start:start + L]


def kmers(sequence, K):
    '''Iterate k-mers for the given sequence.
    '''
    for i in range(len(sequence) - K + 1):
        yield sequence[i:i + K]


