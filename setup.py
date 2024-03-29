from setuptools import setup, find_packages

setup(
    name='pangeapy',
    version='0.0',
    description='initial version',
    
    author='srkim',
    author_email='srkim727@kaist.ac.kr',
    url='https://github.com/srkim727/pangeapy',
    
    packages=find_packages(),
    install_requires=[
        'celltypist', 'scipy'
    ],
)