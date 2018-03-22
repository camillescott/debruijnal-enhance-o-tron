import pytest

from debruijnal_enhance_o_tron.fixtures import ksize, using_ksize

def test_ksize_default(ksize):
    assert ksize == 21


@using_ksize(31)
def test_ksize_override(ksize):
    assert ksize == 31


@using_ksize([25,31])
def test_ksize_override_seq(ksize):
    assert ksize in [25,31]
    if ksize != 25:
        assert ksize == 31
    elif ksize != 31:
        assert ksize == 25
    else:
        assert False
