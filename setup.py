from setuptools import setup

setup(name='qubit_readout',
      version='1.0',
      description='Benchmarking fidelity of single-shot readout in semiconductors',
      author='Daniel Keith',
      author_email='daniel.keith@unsw.edu.au',
      url='https://github.com/silicon-quantum-computing/qubit_readout',
      license = 'Creative Commons for non-commercial use only',
      maintainer          = 'Daniel Keith',
      maintainer_email    = 'daniel.keith@unw.edu.au',
      packages=[ 'qubit_readout', 'qubit_readout.benchmarking'],
     )