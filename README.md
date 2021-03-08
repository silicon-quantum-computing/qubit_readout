# Fidelity

Fidelity is a library for benchmarking the fidelity of single-shot single-spin readout measurements in semiconductors. The overall readout fidelity is used to benchmark the percentage of readout measurements that can correctly identify excited and ground qubit states respectively. Fidelity is intended to comprehensively consider all possible readout errors to set a foundation for a universal and analytical method of benchmarking and comparing readout fidelities across various experiments and systems. The determination of the readout fidelity is broken into two stages:
1. State-to-charge conversion (STC) - where each spin state is uniquely mapped to a directly measurable charge state.
2. Electrical readout (ER) - where each charge state is uniquely identified in the measurement signal.

The fidelity of the STC process can be optimised in terms of the readout time window used to take single-shot measurements. Similarly, the fidelity of the ER process can be optimised in terms of the threshold used to distinguish discrete signal levels that correspond to the qubits charge/spin states. The purpose of the tools within this library is to:
1. Calculate the readout fidelity for a given value of experimental parameters.
2. Optimise the readout fidelity in terms of the readout time and/or signal threshold. 

The tools included here are based on the work presented in ["Benchmarking high fidelity single-shot readout of semiconductor qubits"](https://iopscience.iop.org/article/10.1088/1367-2630/ab242c) that can be refered to for more detailed derivations and discussion.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install the fidelity library from within qubit_readout. Navigate to the directory that contains qubit_readout and run

```bash
pip install setup.py
```

## Usage

Firstly, it is possible to calculate the best readout time to optimise the STC fidelity based on experimental qubit state characteristic transition times. This readout time maximises the probability an excited state transition occurs during the readout time, while minimising the the probability that a ground state transition occurs during the readout time.

```python
import qubit_readout

out_time_excited = 5e-3 #seconds
out_time_ground = 100e-3 #seconds
relax_time = 10 #seconds

qubit_readout.benchmarking.fidelity.optimal_read_time( out_time_excited, out_time_ground, relax_time ) # returns the optimal readout time in seconds 
```

The STC fidelity can be determined in terms of the qubit state characteristic transition times and either the optimised readout time or a given value for the readout time. 

```python
import qubit_readout

out_time_excited = 5e-3 #seconds
out_time_ground = 100e-3 #seconds
relax_time = 10 #seconds

optimal_readout_time = qubit_readout.benchmarking.fidelity.optimal_read_time( out_time_excited, out_time_ground, relax_time )

qubit_readout.benchmarking.fidelity.stc_fidelity( out_time_excited, out_time_ground, relax_time, readout_time = optimal_readout_time ) # returns a STC fidelity for both the ground and excited state
```

Similarly, the ER fidelity can also be determined from the qubit state characteristic transition times, experimental parameters, and a readout time (optimised or not). Instead of finding the optimised threshold for electrical readout the ER fidelity is calculated for a range of normalised thresholds that span the range between two mean signal levels (corresponding to the two possible spin/charge states).
The optimisations for both processes (STC and ER) are done with respect to the readout visibility ($V_{STC/ER}$) which is defined as the sum of the qubit state fidelities ($F^{0}_{STC}$ and $F^{1}_{STC}$, or $F^{0}_{ER}$ and $F^{1}_{ER}$ ), for a given process, minus one i.e.
* $V_{STC} = F^{0}_{STC} + F^{1}_{STC} - 1$.

```python
import qubit_readout

out_time_excited = 5e-3 #seconds
out_time_ground = 100e-3 #seconds
in_time_ground = 10e-3 #seconds
relax_time = 10 #seconds
snr = 5 #signal-to-noise ratio
sample_rate = 10e3 #Hz
filter_cutoff = 5e3 #Hz

optimal_readout_time = qubit_readout.benchmarking.fidelity.optimal_read_time( out_time_excited, out_time_ground, relax_time )

qubit_readout.benchmarking.fidelity.er_fidelity( out_time_excited, out_time_ground, relax_time, readout_time = optimal_readout_time ) 
# returns:
#    - an array of the normalised thresholds the ER fidelity was calculated for,
#    - an array of ER fidelities for the ground state,
#    - an array of ER fidelities for the excited state,
#    - an array of ER visibilities,
#    - the index at which the optimal visbility occurs
```

The overall readout fidelity $F_M$ can then be determined from the STC and ER fidelities as
* $F_M = \frac{F_0 + F_1}{2}$,

where
* $F_0 = F^0_{STC}F^0_{ER} + (1 - F^0_{STC})(1 - F^1_{ER})$, and
* $F_1 = F^1_{STC}F^1_{ER} + (1 - F^1_{STC})(1 - F^0_{ER})$.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)