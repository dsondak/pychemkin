
import os
from setuptools import setup, find_packages

# This reads the __version__ variable from pychemkin/_version.py
exec(open('src/pychemkin/_version.py').read())

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
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    package_data={
        '': [os.path.join('src', 'pychemkin', 'xml_files', '*.xml'),
             os.path.join('src', 'pychemkin', 'thermo_data', '*.sqlite')]
    },
    test_suite='tests'
    )
