"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from os import path
# io.open is needed for projects that support Python 2.7
# It ensures open() defaults to text mode with universal newlines,
# and accepts an argument to specify the text encoding
# Python 3 only projects can skip this import
from io import open

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='arteryfe',  # Required

    version='1.0',  # Required

    description='Implementation of the 1D blood flow equations in FEniCS.',  # Optional

    long_description=long_description,  # Optional
    long_description_content_type='text/markdown',  # Optional (see note above)

    url='https://github.com/akdiem/bloodflow',  # Optional

    author='Syver D. Agdestein, Kristian Valen-Sendstad, Alexandra K. Diem',  # Optional
    author_email='alexandra@simula.no',  # Optional

    classifiers=[  # Optional

        # Pick your license as you wish
        'License :: OSI Approved :: BSD 3-clause',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    packages=find_packages(exclude=['config', 'data', 'docs', 'paper',
                            'tests']),  # Required

    install_requires=['fenics == 2017.2.0', 'numpy', 'h5py'],  # Optional
    extras_require={  # Optional
        'test': ['pytest3'],
    },

)
