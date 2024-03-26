from analysis.mc.interfaces import ImportanceSamplingIntegrator
from typing import Callable, List, Optional

import numpy as np
import pandas as pd
import tensorflow as tf
import tensorflow_probability as tfp
from iflow.integration import integrator
from iflow.integration import couplings
from analysis.nis.masks import binary_masks
from analysis.mc.tools import variance_weighted_result
from analysis.plot_matplotlib import plot_hist
import matplotlib.pyplot as plt

from analysis.import_data import import_true_reco
from analysis.calc import get_kinematics
from typing import Callable

tfd = tfp.distributions
tfb = tfp.bijectors

tf.keras.backend.set_floatx('float64')

def build_nn(in_features, out_features, options, init_debug:bool=False, nhidden:int=4, nchannels:int=256):
    """Build a Dense network with customizable hyperparameters

    Args:
        in_features (_type_): _description_
        out_features (_type_): _description_
        options (_type_): _description_
        init_debug (bool, optional): _description_. Defaults to False.
        nhidden (int, optional): _description_. Defaults to 4.
        nchannels (int, optional): _description_. Defaults to 256.

    Returns:
        _type_: _description_
    """
    del options
    
    bias_initializer = tf.constant_initializer(value=0.01) if init_debug else None
    kernel_initializer = tf.constant_initializer(value=0.01) if init_debug else None

    invals = tf.keras.layers.Input(in_features, dtype=tf.float64)
    hidden = tf.keras.layers.Dense(nchannels, activation='relu', bias_initializer=bias_initializer, kernel_initializer=kernel_initializer)(invals)
    for i in range(nhidden - 1):
        hidden = tf.keras.layers.Dense(nchannels, activation='relu', bias_initializer=bias_initializer, kernel_initializer=kernel_initializer)(hidden)
    #hidden = tf.keras.layers.Dense(32, activation='relu', bias_initializer=bias_initializer, kernel_initializer=kernel_initializer)(hidden)
    #hidden = tf.keras.layers.Dense(32, activation='relu', bias_initializer=bias_initializer, kernel_initializer=kernel_initializer)(hidden)
    outputs = tf.keras.layers.Dense(out_features, bias_initializer='zeros',
                                    kernel_initializer='zeros')(hidden)
    model = tf.keras.models.Model(invals, outputs)
    model.summary()
    return model

def train_iflow(integrate:integrator.Integrator, ptspepoch:int, epochs:int,
                test_callback:Optional[Callable]=None, test_callback_freq:int=50, do_masking:bool=False):
    """ Run the iflow integrator

    Args:
        integrate (Integrator): iflow Integrator class object
        ptspepoch (int): number of points per epoch in training
        epochs (int): number of epochs for training

    Returns:
        numpy.ndarray(float): value of loss (mean) and its uncertainty (standard deviation)

    """
    means = np.zeros(epochs)
    stddevs = np.zeros(epochs)
    losses = np.zeros(epochs)
    for epoch in range(epochs):
        loss, integral, error = integrate.train_one_step(ptspepoch, integral=True, do_masking=do_masking)
        means[epoch] = integral
        stddevs[epoch] = error
        losses[epoch] = loss
        _, current_precision = variance_weighted_result(means[:epoch+1], stddevs[:epoch+1])
        if epoch % 10 == 0:
            print('Epoch: {:3d} Loss = {:8e} Integral = '
                  '{:8e} +/- {:8e} Total uncertainty = {:8e}'.format(epoch, loss,
                                                                     integral, error,
                                                                     current_precision))
        
        if epoch % test_callback_freq == 0 and test_callback is not None:
            test_callback(integrate)

    return means, stddevs, losses

def build_flow_dist(boundaries:np.ndarray, nbins:int=16, build_nn:Callable=build_nn):
    # Calculate scale and hist
    ndims = len(boundaries)
    width = (boundaries[:, 1] - boundaries[:, 0])
    shift = (boundaries[:, 1] + boundaries[:, 0])/2 - width/2

    base = tfd.Uniform(low=np.zeros(ndims, dtype=np.float64), high=np.ones(ndims, dtype=np.float64))
    bijectors = []

    for mask in binary_masks(ndims):
        bijectors.append(couplings.PiecewiseRationalQuadratic(mask, build_nn, num_bins=nbins, options=None))

    bijectors.append(tfp.bijectors.Scale(width))
    bijectors.append(tfp.bijectors.Shift(shift))

    bijector = tfb.Chain(list(reversed(bijectors)))

    return tfd.TransformedDistribution(distribution=tfd.Independent(distribution=base, reinterpreted_batch_ndims=1), bijector=bijector)

def plot_marginals(dist, n_samples:int=10000, is_dist:bool=False, separate:bool=True, plot_args:dict={}, plot_grid=(2,4), plot_size=(20, 10)):
    samples = dist.sample(n_samples).numpy()
    plt_args={ 'scientific_stats': True, **plot_args }
    
    if not separate:
        all_max = np.max(samples).max()
        all_min = np.min(samples).min()
        dfbase = {}
        for i in range(len(samples.T)):
            dfbase[f'x{i}'] = samples.T[i]

        return plot_hist(pd.DataFrame(dfbase), normalize=True, xlim=(all_min, all_max), xlim_binning=(all_min, all_max), **plt_args)
    else:
        fig, axes = plt.subplots(nrows=plot_grid[0], ncols=plot_grid[1], figsize=plot_size)
        for i in range(samples.shape[1]):
            plot_hist({ f'Dim {i+1}': samples.T[i] }, normalize=True, ax=axes.flat[i], **plt_args)
        
        to_del = int(np.prod(plot_grid) - len(samples.T))
        for j in range(to_del):
            fig.delaxes(axes.flat[-(j+1)])
        
        return fig, axes

def build_mem_integrand(reco_kin, print_ps_efficiciency:bool=True)->Callable:
    from analysis.cffi.mg5.lib import mc_batch, mc_batch_sigma
    
    def integrand(vars:np.ndarray, mode=1, me_type=1):
        #return np.sum(vars, axis=-1)
        
        if print_ps_efficiciency:
            tf.print(f"PS points given:found [{len(vars)}:", end="")
        
        res = mc_batch(reco_kin, vars.flatten().tolist(), mode=mode, me_type=me_type)
        
        found = np.count_nonzero(res)
        efficiency = found/len(vars)
        
        if print_ps_efficiciency:
            tf.print(f"{found}] ({(100*efficiency):6.2f} %) ")
        
        return tf.convert_to_tensor(res, dtype=tf.float64), tf.convert_to_tensor(efficiency, dtype=tf.float64)
    
    return integrand

class IFlowIntegrator(ImportanceSamplingIntegrator):
    def __init__(self,
                 integrand:Callable,
                 boundaries:List[List[float]],
                 optimizer_args:dict={}, integrator_args:dict={},
                 *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
    
        optim_args = { 'clipnorm': 10.0,  **optimizer_args }
        optimizer = tf.keras.optimizers.Adam(lr, **optim_args)
        
        integr_args = { 'loss_func': 'exponential', **integrator_args }
        integrate = integrator.Integrator(func, dist, optimizer, **integr_args)