import pytest
from collections import OrderedDict
import inspect

from debruijnal_enhance_o_tron.sequence import kmers


class Collector(object):

    def __init__(self):
        self.items = OrderedDict()
    
    def add(self, src_name, item):
        if src_name in self.items:
            self.items[src_name].append(item)
        else:
            self.items[src_name] = [item]

    def pop(self):
        self.items.popitem()

    def clear(self):
        self.items.clear()

    def values(self):
        for key, item_list in self.items.items():
            for item in item_list:
                yield item


@pytest.fixture
def consume_collector():
    seq_collector = Collector()

    def add(*args):
        for sequence in args:
            seq_collector.add(inspect.stack()[1][3], sequence)
        return seq_collector

    return add


@pytest.fixture
def consume(request, consume_collector, graph):
    seq_collector = consume_collector()
   
    def _consume():
        print(seq_collector.items)
        for sequence in seq_collector.values():
            graph.add(sequence)
    
    return _consume

@pytest.fixture
def check_fp_collector(request):
    check_collector = Collector()

    def add(*args):
        for func_pair in args:
            if type(func_pair) is not tuple:
                raise TypeError('_check_fp_collector expects a tuple of (function, result)')
            check_collector.add(inspect.stack()[1][3], func_pair)
        return check_collector

    return add


def check_fp_xfail(request):
    pytest.xfail('False-positive check FAILED ({}).'.format(request.fixturename))


def check_fp_pass(request):
    print('False-positive check PASS ({}).'.format(request.fixturename))


@pytest.fixture
def check_fp(request, check_fp_collector, consume_collector, graph, ksize):
    
    def _check_fp():
        _graph = graph.shallow_clone()
        funcs = check_fp_collector()
        sequences = consume_collector()

        for sequence in sequences.values():
            for kmer in kmers(sequence, ksize):
                if _graph.get(kmer):
                    check_fp_xfail(request)
        
        for sequence in sequences.values():
            _graph.add(sequence)

        results = [func(_graph) == result for func, result in funcs.values()]
        if not all(results):
            check_fp_xfail(request)
        else:
            check_fp_pass(request)
    
    return _check_fp
