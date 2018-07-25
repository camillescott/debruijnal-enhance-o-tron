import pytest

from debruijnal_enhance_o_tron.sequence import (mutate_base,
                                                mutate_sequence,
                                                mutate_position,
                                                get_random_sequence,
                                                reads,
                                                kmers,
                                                revcomp)

@pytest.fixture(params=[21, 51], ids=lambda k: 'K={0}'.format(k))
def ksize(request):
    '''
    Core fixture to set the K parameter.
    '''
    return request.param


def using_ksize(K):
    '''
    Convenience wrapper to reduce verbosity of indirect parametrization
    '''
    def wrapped(fixture_func):
        if isinstance(K, int):
            ksize = [K]
        else:
            ksize = list(K)
        return pytest.mark.parametrize('ksize', 
                                       ksize,
                                       indirect=['ksize'],
                                       ids=lambda k: 'K={0}'.format(k))(fixture_func)
    return wrapped


@pytest.fixture(params=[500, 1000], ids=lambda l: 'L={0}'.format(l))
def length(request):
    '''
    Basic length fixture.
    '''
    return request.param


def using_length(L):
    '''
    Convenience wrapper for basic length fixure.
    '''

    def wrapped(fixture_func):
        if isinstance(L, int):
            length = [L]
        else:
            length = list(L)
        return pytest.mark.parametrize('length', 
                                       length,
                                       indirect=['length'],
                                       ids=lambda l: 'L={0}'.format(l))(fixture_func)
    return wrapped


@pytest.fixture
def random_sequence(request, ksize, length):
    global_seen = set()

    def get(exclude=None, length=length):

        try:
            sequence = get_random_sequence(length, 
                                           ksize,
                                           exclude=exclude,
                                           seen=global_seen)
        except ValueError:
            pytest.xfail('Failed to generate sequence of given length and k-size')

        for kmer in kmers(sequence, ksize-1):
            global_seen.add(kmer)

        return sequence

    return get

