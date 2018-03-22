from functools import wraps

import pytest

@pytest.fixture(params=[21], ids=['K=21'])
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
