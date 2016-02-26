from setuptools import setup

setup(name='experimentdataanalysis',
    version='0.1',
    description='parsing data output from experiment and performing analysis',
    url='http://github.com/mwmacmahon/experimentdataanalysis',
    author='Michael Macmahon',
    author_email='mwmacmahon@gmail.com',
    setup_requires=['pytest-runner'],
    tests_require=['pytest']
)
