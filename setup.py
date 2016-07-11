# With reference to
# https://jeffknupp.com/blog/2013/08/16/open-sourcing-a-python-project-the-right-way/

from __future__ import print_function
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import io
import os
import sys

import experimentdataanalysis

here = os.path.abspath(os.path.dirname(__file__))


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.txt', 'CHANGES.txt')


class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

#    def finalize_options(self):
#        TestCommand.finalize_options(self)
#        self.pytest_args = []
#        self.test_suite = True

    def run_tests(self):
        import pytest
        if __name__ == "__main__":
            errcode = pytest.main(self.pytest_args)
            sys.exit(errcode)

setup(
    name='experimentdataanalysis',
    version=experimentdataanalysis.__version__,
    description='parsing data output from experiment and performing analysis',
    long_description=long_description,
    url='http://github.com/mwmacmahon/experimentdataanalysis',
    author='Michael Macmahon',
    author_email='mwmacmahon@gmail.com',
    tests_require=['pytest'],     # seemed to have fixed multiprocessing issues
    cmdclass = {'test': PyTest},  # but can always run as python -m py.test
    packages=find_packages(),
    include_package_data=True,

#   TODO: GIVES AN ERROR, FIX:
#    install_requires=[i.strip() for i in open("requirements.txt").readlines()],

#    license='',
#    setup_requires=[],
#    platforms='any',
#    test_suite='experimentdataanalysis.test.test_experimentdataanalysis',
#    scripts = [],
#    extras_require={
#        'testing': ['pytest'],
#    }

    #  # ONLY FOR SCRIPTS INSIDE MODULES, NOT FOR SCRIPT-OUTSIDE-OF-MODULE 
    #entry_points={
    #    'console_scripts': [
    #        'experimentdataanalysis_databrowser = testscript_databrowser:main',
    #        'experimentdataanalysis_process1dtrkr = testscript_process1Dtrkr:main',
    #        'experimentdataanalysis_process3dtrkr = testscript_process3Dtrkr:main',
    #    ],
    #    'gui_scripts': [
    #        'experimentdataanalysis_databrowser_gui = testscript_databrowser:main',
    #    ]
    #}
)