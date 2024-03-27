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

DEFAULT_HIDDEN_LAYERS = 4
DEFAULT_HIDDEN_CHANNELS = 256

def build_nn(in_features, out_features, options, init_debug:bool=False, nlayers:Optional[int]=None, nchannels:Optional[int]=None):
    """Build a Dense network with customizable hyperparameters

    Args:
        in_features (_type_): _description_
        out_features (_type_): _description_
        options (_type_): _description_
        init_debug (bool, optional): _description_. Defaults to False.
        nlayers (int, optional): _description_. Defaults to 4.
        nchannels (int, optional): _description_. Defaults to 256.

    Returns:
        _type_: _description_
    """
    del options
    
    nlayers = DEFAULT_HIDDEN_LAYERS if nlayers is None else nlayers
    nchannels = DEFAULT_HIDDEN_CHANNELS if nchannels is None else nchannels
    
    bias_initializer = tf.constant_initializer(value=0.01) if init_debug else None
    kernel_initializer = tf.constant_initializer(value=0.01) if init_debug else None

    invals = tf.keras.layers.Input(in_features, dtype=tf.float64)
    hidden = tf.keras.layers.Dense(nchannels, activation='relu', bias_initializer=bias_initializer, kernel_initializer=kernel_initializer)(invals)
    for i in range(nlayers - 1):
        hidden = tf.keras.layers.Dense(nchannels, activation='relu', bias_initializer=bias_initializer, kernel_initializer=kernel_initializer)(hidden)
        
    #hidden = tf.keras.layers.Dense(32, activation='relu', bias_initializer=bias_initializer, kernel_initializer=kernel_initializer)(hidden)
    #hidden = tf.keras.layers.Dense(32, activation='relu', bias_initializer=bias_initializer, kernel_initializer=kernel_initializer)(hidden)
    outputs = tf.keras.layers.Dense(out_features, bias_initializer='zeros',
                                    kernel_initializer='zeros')(hidden)
    model = tf.keras.models.Model(invals, outputs)
    model.summary()
    return model

def build_flow_dist(boundaries:np.ndarray, nbins:int=16, build_nn:Callable=build_nn, nn_args:dict={}):
    # Calculate scale and hist
    ndims = len(boundaries)
    width = (boundaries[:, 1] - boundaries[:, 0])
    shift = (boundaries[:, 1] + boundaries[:, 0])/2 - width/2

    base = tfd.Uniform(low=np.zeros(ndims, dtype=np.float64), high=np.ones(ndims, dtype=np.float64))
    bijectors = []

    nn_args_in = {  }

    for mask in binary_masks(ndims):
        bijectors.append(couplings.PiecewiseRationalQuadratic(mask, build_nn, num_bins=nbins, options=None))

    bijectors.append(tfp.bijectors.Scale(width))
    bijectors.append(tfp.bijectors.Shift(shift))

    bijector = tfb.Chain(list(reversed(bijectors)))

    return tfd.TransformedDistribution(distribution=tfd.Independent(distribution=base, reinterpreted_batch_ndims=1), bijector=bijector)

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

def plot_marginals(dist, n_samples:int=10000, separate:bool=True, plot_args:dict={}, plot_grid=(2,4), plot_size=(20, 10)):
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
        
        return fig

def plot_summary(losses:np.ndarray, means:np.ndarray, stddevs:np.ndarray):
    from analysis.plot_matplotlib import get_colorpalette

    yscale = 'log'

    palette = get_colorpalette()
    fig, ax = plt.subplots(3, 1, figsize=(4,8), gridspec_kw={'hspace': 0, 'wspace': 0})
    ax[0].plot(np.arange(len(losses)), losses, label=r'Loss', color=palette[0])
    ax[0].set_xlabel(r'Epoch')

    ax[1].plot(np.arange(len(means)), means, label=r'$I$', color=palette[1])
    ax[2].plot(np.arange(len(stddevs)), stddevs,label=r'$\sigma$', color=palette[2])

    for i in range(3):
        ax[i].set_yscale(yscale)

    for i in range(len(ax)):
        ax[i].legend(loc='best')
        
    return fig

class IFlowIntegrator(ImportanceSamplingIntegrator):
    def __init__(self,
                 integrand:Callable,
                 boundaries:List[List[float]],
                 lr:float=1e-3, optimizer_args:dict={}, integrator_args:dict={},
                 nbins:int=30, build_nn:Callable=build_nn,
                 nlayers:int=4, nchannels:int=256,
                 *args, **kwargs) -> None:
        super().__init__(integrand=integrand, boundaries=boundaries)
        
        DEFAULT_HIDDEN_CHANNELS = nchannels
        DEFAULT_HIDDEN_LAYERS = nlayers
        
        dist = build_flow_dist(np.array(boundaries), nbins=nbins, build_nn=build_nn)
    
        optim_args = { 'clipnorm': 10.0,  **optimizer_args }
        optimizer = tf.keras.optimizers.Adam(lr, **optim_args)
        
        integr_args = { 'loss_func': 'exponential', **integrator_args }
        integrate = integrator.Integrator(integrand, dist, optimizer, **integr_args)
        
        self.dist = dist
        self.iflow = integrate
        self.optimizer = optimizer
        
        self.means = []
        self.stddevs = []
        self.losses = []
        
    def sample(self, nsamples:int):
        return self.iflow.sample(nsamples)        
    
    @tf.function(input_signature=[tf.TensorSpec(None, tf.double)])
    def call_func(self, vals):
        """Allows calling external functions"""
        return tf.numpy_function(self.integrand, [vals], [tf.double, tf.double])
        
    @tf.function
    def step(self, nsamples, integral:bool=False, do_masking:bool=False):
        samples = self.dist.sample(nsamples)
        # self.samples = tf.concat([self.samples, samples], 0)
        # if self.samples.shape[0] > 5001:
        #     self.samples = self.samples[nsamples:]
        results = self.call_func(samples)
        true, efficiency = tf.abs(results[0]), results[1]
        true.set_shape([None])
        #tf.print('Efficiency', efficiency)
        
        with tf.GradientTape() as tape:
            test = self.dist.prob(samples)
            logq = self.dist.log_prob(samples)
            mean, var = tf.nn.moments(x=true/test, axes=[0])
            true = tf.stop_gradient(true/mean)
            logp = tf.where(true > 1e-16, tf.math.log(true),
                            tf.math.log(true+1e-16))
            
            # mask zero PS items
            if do_masking:
                mask = tf.cast(true > 0, dtype=tf.bool)  
                mask.set_shape([None])

                true = tf.boolean_mask(true, mask)
                test = tf.boolean_mask(test, mask)
                logp = tf.boolean_mask(logp, mask)
                logq = tf.boolean_mask(logq, mask)
            
            # loss = self.loss_func(samples, samples, 1e-1, true, test, 100)
            #loss = tf.stop_gradient(-tf.math.log(efficiency))*self.loss_func(true, test, logp, logq)
            loss = self.iflow.loss_func(true, test, logp, logq)

        grads = tape.gradient(loss, self.iflow.dist.trainable_variables)
        self.optimizer.apply_gradients(
            zip(grads, self.iflow.dist.trainable_variables))

        if integral:
            return loss, mean, tf.sqrt(var/(nsamples-1.))

        return loss
    
    def train(self, ptspepoch:int, epochs:int, test_callback:Optional[Callable]=None, test_callback_freq:int=50, do_masking:bool=False):
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
            loss, integral, error = self.step(ptspepoch, integral=True, do_masking=do_masking)
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
                test_callback(self)

        self.means = np.concatenate([self.means, means])
        self.stddevs = np.concatenate([self.stddevs, stddevs])
        self.losses = np.concatenate([self.losses, losses])

        return means, stddevs, losses
    
    def integrate(self, sample_size:int=5000, nepochs:int=500, do_masking:bool=False):
        self.train(ptspepoch=sample_size, epochs=nepochs, do_masking=do_masking)
        
        tot_res, tot_uncert = variance_weighted_result(self.means, self.stddevs)
        
        return tot_res, tot_uncert
    
    def plot_marginals(self, **kwargs):
        return plot_marginals(self.iflow.dist, **kwargs)
    
    def plot_summary(self):
        return plot_summary(self.losses, self.means, self.stddevs)