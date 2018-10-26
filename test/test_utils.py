import sys
sys.path.insert(0, 'arteryfe/')

from utils import *
from utils import is_near as near


def test_unit_to_mmHg():
    for p in range(0, 30000, 500):
        print(unit_to_mmHg(p), p/1333.22387415)
        assert(near(unit_to_mmHg(p), p/1333.22368421053))


def test_mmHg_to_unit():
    for p in range(0, 300, 5):
        print(mmHg_to_unit(p), p*1333.22387415)
        assert(near(mmHg_to_unit(p), p*1333.223684210535))


def test_adimensionalise():
    None


def test_redimensionalise():
    None


def test_read_output():
    None


def test_XDMF_to_matrix():
    None


def test_plot_matrix():
    None
