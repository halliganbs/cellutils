from importlib.metadata import entry_points
from setuptools import setup

import os

if os.name == 'nt':
    os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='cellutils',
    version='0.0.2',
    description='High throughput microscopy functions',
    long_description=readme(),
    classifiers=[
        'Development Status :: Number - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.10'
    ],
    keywords='mircoscopy',
    url="https://github.com/halliganbs/cellutils",
    author='halliganbs',
    author_email='bhalliga@med.umich.edu',
    license='MIT',
    packages=['cellutils']
)