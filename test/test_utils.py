import sys

sys.path.insert(0, 'src/')
from utils import *
from utils import is_close as near


def test_unit_to_mmHg():
	for p in range(0, 30000, 500):
		assert(near(unit_to_mmHg(p), p/1333.22387415))


def test_mmHg_to_unit():
	for p in range(0, 300, 5):
		assert(near(mmHg_to_unit(p), p/1333.22387415))


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


def test_utils(config_location)

	test_unit_to_mmHg():

	test_mmHg_to_unit():

	test_adimensionalise():

	test_redimensionalise():

	test_read_output():

	test_XDMF_to_matrix():

	test_plot_matrix():


if __name__ == '__main__'
	test_utils(sys.argv[1])
