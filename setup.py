
import os
from setuptools import setup, find_packages

# This reads the __version__ variable from pychemkin/_version.py
exec(open('pychemkin/_version.py').read())

# Get the long description from the README file
long_description = open('README.md').read()

# Read in requirements.txt
requirements = open('requirements.txt').readlines()
requirements = [r.strip() for r in requirements]

setup(
    name='pychemkin',
    version=__version__,
    description='Chemical kinetics library',
    keywords='chemical kinetics',
    long_description=long_description,
    install_requires=requirements,
    url='https://github.com/dsondak/pychemkin',
    author='PyChemKin Developers',
    author_email='hsim13372@gmail.com',
    license='MIT',
    packages=['pychemkin'],
    test_suite='tests',
    include_package_data=True,
    package_data={
        'pychemkin': [os.path.join('pychemkin', 'thermo_database', '*.sqlite')]
    }
    )
