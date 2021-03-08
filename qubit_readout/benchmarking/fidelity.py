# Copyright (C) 2021 by “Silicon Quantum Computing Pty. Ltd.”
# This software/code is provided by the copyright owner “Silicon Quantum Computing” and may be distributed
# under the terms of the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public
# License ("Public License") – https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode.
# By exercising the Licensed Rights (defined in the link provided above), you accept and agree to be bound by the terms
# and conditions of this Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License ("Public License").

# A license for commercial purposes may be available and obtained by contacting the copyright holder at info@sqc.com.au.

# NO WARRANTY
# THIS SOFTWARE/CODE IS PROVIDED BY THE COPYRIGHT HOLDER "SILICON QUANTUM COMPUTING" AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF NONINFRINGEMENT, MERCHANTIBILITY
# AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
# THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR SPECIAL, EXEMPLARY,
# OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGES.

# Author - Daniel Keith <daniel.keith@unsw.edu.au> 2021

import functools
from functools import lru_cache
import math
from numba import jit, float32
import numpy
import operator
import scipy
from scipy import signal

def optimal_read_time( out_time_excited, out_time_ground, relax_time ):
    '''
    Calculate optimal readout time for single-shot spin readout given characteristic transition times.

    Parameters
    ----------
    out_time_excited : float
        Characteristic transitioin time in seconds of the EXCITED spin state OUT of the quantum dot to the reservoir.
    out_time_ground : float
        Characteristic transition time in seconds of the GROUND spin state OUT of the quantum dot to the reservoir.
    relax_time : float
        Relaxation time in seconds of the excited spin state to the ground state.

    Returns
    -------
    read_time : float
        Optimal time window in seconds for single-shot readout.
    '''

    x = relax_time * ( out_time_ground - out_time_excited ) + out_time_excited * out_time_ground
    read_time = ( relax_time * out_time_ground * out_time_excited / x ) * math.log( out_time_ground * ( relax_time + out_time_excited ) / ( relax_time * out_time_excited ) )
    return read_time

def stc_fidelity( out_time_excited, out_time_ground, relax_time, readout_time = None ):
    '''
    Calculate optimal spin state fidelities for single-shot spin readout for the given state transition times and readout time.

    Parameters
    ----------
    out_time_excited : float
        Characteristic transition time in seconds of the EXCITED spin state OUT of the quantum dot to the reservoir.
    out_time_ground : float 
        Characteristic transition time in seconds of the GROUND spin state OUT of the quantum dot to the reservoir.
    relax_time : float
        Relaxation time in seconds of the excited spin state to he ground state.
    readout_time : float, optional 
        Readout time window in seconds in which a transition event can be detected. Default is None which will automatically optimise the readout time based on the provided characteristic transition times and relaxation time.
    
    Returns
    -------
    ground_fid : float 
        Probability that a transition correctly does not occur during each individual readout time window when the qubit is initially in the ground state.
    excited_fid : float 
        Probability that a transition correctly does occur during each indivdual readout time window when the qubit is initially in the excited state.
    
    '''

    if readout_time is None: readout_time = optimal_read_time( out_time_excited, out_time_ground, relax_time )

    x = relax_time * ( out_time_ground - out_time_excited ) + out_time_excited * out_time_ground
    excited_fid = ( 1 / x ) * ( ( 1 - math.exp( -1*readout_time / out_time_ground ) ) * out_time_ground * out_time_excited + 
            ( math.exp( - ( ( relax_time + out_time_excited ) * readout_time ) / ( out_time_excited * relax_time )) - 1 ) * relax_time * ( out_time_excited - out_time_ground ) )
    ground_fid = math.exp( -1*readout_time / out_time_ground )
    
    return ground_fid, excited_fid

def er_fidelity( out_time_excited, in_time_ground, readout_time, snr, sample_rate, filter_cutoff, threshold_num=1001 ):
    '''
    Calculates electrical fidelities for the given state characteristic transition times and experimental parameters for a range of threshold values and returns the optimal threshold and correspopnding electrical fidelities.

    Parameters
    ----------
    out_time_excited : float
        Characteristic transition time in seconds of the EXCITED spin state OUT of the quantum dot to the reservoir.
    in_time_ground : float
        Characteristic transition time in seconds of the GROUND spin state IN to the quantum dot from the reservoir.
    readout_time : float
        Duration in seconds of each single-shot measurement used to detect state transition events.
    snr : float
        Average voltage signal-to-noise ratio between the single-shot readout measurements corresponding to the charge states.
    sample_rate : float
        Data points taken per second.
    filter_cutoff : 
        Frequency in Hertz corresponding to the -3dB point of the limiting filter used during single-shot readout measurements.
    threshold_num : int, optional
        Number of thresholds to test for the optimal electrical fidelity. Thresholds are linearly spaced over a normalised range corresponding to the values between the higher and lower signal levels measured.

    Returns
    -------
    thresholds : 1d array
        Thresholds used to calculate that the electrical readout fidelity.
    ground_fid : 1d array
        Probabilities of no transition being detected within the readout time window when no transition occurred for each threshold value tested. 
    excited_fid : 1d array
        Probabilities of a transition being detected within the readout time window when a transition occurred for each threshold value tested. 
    vis : 1d array
        Visibility defined as ground_fid + excited_fid - 1 for each threshold value tested. The visibility is used to optimise the electrical fidelities.
    idx : int 
        The index of thresholds corresponding to the threshold value in which the optimal visibility, and hence fidelities, occurs. 
    '''

    @lru_cache
    @jit(float32(float32, float32, float32))
    def ncdf(x, mu, sigma):
        r_val = 0.5 * (1 + math.erf((x - mu) / (numpy.sqrt(2) * sigma)))
        return r_val
    ncdf2 = numpy.vectorize(ncdf)

    sample_time = 1 / sample_rate

    fr = 2 * filter_cutoff * sample_time
    tnr = 2 * fr / (fr + 1) if fr < 1 else 1
    
    nr = ( readout_time / sample_time ) * tnr  
    nh = ( in_time_ground / sample_time )  
    nl = ( out_time_excited / sample_time )  

    mul = 0
    muh = 1

    noise = 1 / snr

    thresholds = numpy.linspace(0,1,threshold_num) 

    cl = numpy.power((ncdf2(thresholds, mul, noise)), nr)  

    b, a = scipy.signal.bessel(8, 2 * numpy.pi * filter_cutoff, 'low', analog=True, norm='phase')

    @lru_cache
    def funa(f):
        num = numpy.polyval(b, f * 1j)
        den = numpy.polyval(a, f * 1j)
        return numpy.abs(num / den)

    nr_inv = 1 / nr
    t_inv = 1 / ( tnr * sample_time )
    ncdf_l = ncdf2( thresholds, mul, noise )

    nl_inv = 1 / nl
    nh_inv = 1 / nh

    @lru_cache
    def Sn(n, v, vj):
        vj = int(vj)
        part1 = n * nr_inv * (ncdf2(v, (muh - mul) * funa(t_inv / n) + mul, noise))
        part2 = (1 - n * nr_inv) * ncdf_l[vj]
        r_val = numpy.power((part1 + part2), nr)
        return r_val

    @lru_cache
    @jit(float32(float32))
    def en(s):
        return (numpy.exp((1 - s) * nl_inv)) / (1 - numpy.exp(- nr * nl_inv)) * nl_inv

    @lru_cache
    @jit(float32(float32))
    def etan(n):
        return (numpy.exp((1 - n) * nh_inv)) * nh_inv

    @lru_cache
    def fun1(v, j, n, s):
        return en(s) * etan(n) * Sn(n, v, j)

    @lru_cache
    def fun2(v, j, n, s):
        return en(s) * etan(n) * Sn(nr - s, v, j)

    ch = numpy.array( [ scipy.integrate.dblquad( functools.partial( fun1, v_val, j ), 1, nr - 1, lambda x: 1, lambda x: nr - x)[0] + 
                        scipy.integrate.dblquad( functools.partial( fun2, v_val, j ), 1, nr - 1, lambda x: nr - x, lambda x: 20 * nh)[0] for j, v_val in enumerate(thresholds) ] )

    ri = tnr / nh  
    ro = tnr / nl  
    if ro == ri:
        ro += 0.0001

    navg = 1 - (((1 - numpy.exp(ro / 2 - ri / 2)) * ro) / ((1 - numpy.exp(ro / 2)) * (ro - ri)))

    ground_fid = cl
    excited_fid = (1 - navg) * (1 - ch) + navg * (1 - cl)

    vis = ground_fid + excited_fid - 1

    idx = numpy.argmax(vis)

    return thresholds, ground_fid, excited_fid, vis, idx        