import gc
from torch import nn
import torch
import numpy as np

from analysis.utils import module_reload
module_reload('normflows')
import normflows as nf

from analysis.mc import ImportanceSamplingIntegrator
from analysis.nis.distributions import HyperUniform, Uniform
from analysis.nis.masks import binary_masks
from analysis.nis.transforms import ScaleShift
from analysis.nis.divergences import kl_div, exp_div
from analysis.nis.couplings import PiecewiseRationalQuadraticSpline
from analysis.nis.splines import rational_quadratic_spline, DEFAULT_MIN_BIN_HEIGHT, DEFAULT_MIN_BIN_WIDTH, DEFAULT_MIN_DERIVATIVE

from typing import Optional, Union, Literal, Callable, Union
from tqdm.auto import tqdm

def build_param_net(hidden_channels:int, num_bins:int, n_layers:int=5,
                    activation:Callable=torch.nn.ReLU, dropout_probability:float=0., debug:bool=False,
                    bay_nsf:bool=False):
    
    def builder(in_features , out_features, options:None):
        """in_features may be features to condition on

        Args:
            in_features (_type_): _description_
            out_features (_type_): _description_

        Returns:
            _type_: _description_
        """
        class DenseNetwork(nn.Sequential):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                
                hidden_layers = []
                for i in range(n_layers-1):
                    hidden_layers.append(nn.Linear(hidden_channels if i > 0 else in_features, hidden_channels, bias=True))
                    hidden_layers.append(activation())
                    
                self.hidden_layers = nn.ModuleList(hidden_layers)
                self.final_layer = nn.Linear(hidden_channels if n_layers > 1 else in_features, out_features, bias=True)
            
            def forward(self, z, context=None):
                for layer in self.hidden_layers:
                    z = layer(z)
                
                return self.final_layer(z)
        
        seq = DenseNetwork()
        
        if debug:
            for hl in seq.hidden_layers:
                if isinstance(hl, nn.Linear):
                    nn.init.constant_(hl.weight, 0.01)
                    nn.init.constant_(hl.bias, 0.01)
            
        nn.init.zeros_(seq.final_layer.weight)
        nn.init.zeros_(seq.final_layer.bias)
        
        if bay_nsf:
            nn.init.constant_(seq.final_layer.bias[2*num_bins:], np.log(np.exp(1 - DEFAULT_MIN_DERIVATIVE) - 1))
        
        return seq
    
    return builder
    

class NeuralImportanceSamplingIntegrator(ImportanceSamplingIntegrator):
    def __init__(self,
                 integrand:Callable,
                 boundaries:list[list[float]],
                 integrand_format:Literal['numpy','list','torch']='torch',
                 
                 hidden_channels:int=256,
                 hidden_layers:int=2,

                 # Rational Spline Coupling options
                 num_bins:int=16,
                 hypercube_width:float=1.,
                 init_identity:bool=True,
                 seed:Optional[int]=43,
                 cuda_if_available:bool=True,
                 device:Optional[Union[torch.device,str]]=None,
                 optimizer_args:dict={},
                 lr:float=1e-3,
                 silent:bool=False,
                 shift_scale:bool=True,
                 debug:bool=False,
                 *args, **kwargs):
        
        super().__init__(integrand=integrand, boundaries=boundaries)
        
        device = (torch.device(device) if isinstance(device, str) else device) if device is not None else (
            torch.device('cuda' if torch.cuda.is_available() and cuda_if_available else 'cpu'))
        
        boundaries = np.array(self.boundaries)
        input_size = len(boundaries.T)
            
        # Ensure couplings allow information flow through feature dimensions
        masks = binary_masks(input_size)

        if seed is not None:
            torch.manual_seed(seed)
            
        flow_args = { 'num_bins': num_bins, 'tails': None, #'tails': 'linear', 'tail_bound': hypercube_width/2,
                     'init_identity': init_identity, 'device': device, }
                     #'build_transform_net_create_fn': build_param_builder }
        flows = []
        i = 0
        for mask in masks:
            #flows += [nf.flows.CoupledRationalQuadraticSpline(input_size, hidden_layers, hidden_units,
            #                                                  reverse_mask=not (not 1 % 2), feature_mask=torch.tensor(mask, dtype=torch.int32), **flow_args )]
            flows += [PiecewiseRationalQuadraticSpline(mask, build_param_net(hidden_channels=hidden_channels, num_bins=num_bins, n_layers=hidden_layers, debug=debug),
                                                       num_bins=num_bins, i_test=i )]
            i += 1

        # Set base distribuiton
        if shift_scale:            
            scale = torch.tensor((boundaries[1] - boundaries[0]), dtype=torch.double)
            shift = torch.tensor((boundaries[1] + boundaries[0])/2 - (boundaries[1] - boundaries[0])/2, dtype=torch.double)
            
            flows += [ScaleShift(scale=scale, shift=shift)]
            
            
        # Construct flow model
        #self.q0 = Uniform(input_size, low=-hypercube_width/2, high=hypercube_width/2, device=device)
        optim_args = { 'lr': lr, } #'weight_decay': 1e-3 }
        self.optim_args = {**optim_args, **optimizer_args}
        
        self.n_couplings = len(masks)  
        self.q0 = Uniform(input_size, low=0., high=1., device=device, deterministic=debug)
        self.nfm = nf.NormalizingFlow(q0=self.q0, flows=flows).to(device)
        self.optimizer = torch.optim.Adam(self.nfm.parameters(), **self.optim_args)
        self.call_scheduler = lambda a: None
        self.integrand_format = integrand_format
        self.silent = silent
        self.device = device
        
        self.best_params = None
        self.best_epoch = None
        self.best_sigma = 1e40
        self.last_n_worse = 0
        self.n_resets = 0
        self.using_best = False
        
        self.epoch_abs = 0
        self.statistics = {
            'train_results': [],
            'train_losses': [],
            'train_sigmas': []
        }
    
    def step(self, loss_func:Optional[Callable]=kl_div, n_samples:int=64000, with_loss:bool=True, update_best:bool=True):            
        samples, logq = self.nfm.sample(n_samples)
        samples = samples.detach()
        
        # Move to CPU and convert to list/numpy array, if necessary
        if self.integrand_format == 'numpy':
            samples = samples.cpu().numpy()
            true = np.abs(self.integrand(samples), dtype=np.double)
        elif self.integrand_format == 'list':
            samples = samples.tolist()
                
        true = torch.as_tensor(self.integrand(samples), dtype=torch.double, device=self.device).abs()
        #with torch.no_grad():
        test = logq.exp()
        
        mean = (true/test).detach().mean()
        var = (true/test).detach().var().item()
        
        true = true/mean
     
        logp = torch.where(true > 1e-16, true.log(), (true + 1e-16).log())
        #print('avg_logp', logp.mean())
        
        sigma = np.sqrt(var/(n_samples-1))
        
        print(f'true {(true.sum()).item():.3f}')
        print(true.shape, test.shape, logp.shape, logq.shape)
        
        #self.epoch_abs += 1
        self.statistics['train_results'].append(mean.item())
        self.statistics['train_sigmas'].append(sigma)
        
        if update_best:
            if sigma < self.best_sigma:
                self.best_sigma = sigma
                self.best_params = self.nfm.state_dict()
                self.best_epoch = self.epoch_abs
                self.last_n_worse = 0
                self.using_best = True
            else:
                self.using_best = False
                self.last_n_worse += 1
        
        if with_loss:
            loss = loss_func(true, test, logp, logq)
            
            self.statistics['train_losses'].append(loss.item())
            return loss, mean, var, sigma
        else:
            self.statistics['train_losses'].append(None)
            return mean, var, sigma
        
    def adapt(self, loss_func:Callable=exp_div, n_samples:int=64000,
              target_prec:float=1e-3, n_epochs_max:int=10, test_callback:Optional[Callable]=None,
              n_max_fail:Optional[int]=3, ):
        
        tot_precision = 1e90 if self.epoch_abs == 0 else variance_weighted_result(np.array(self.statistics['train_results']),
                                                                                        np.array(self.statistics['train_sigmas']))[1]
        
        self.nfm.train()
        for epoch in (pbar := tqdm(range(n_epochs_max), disable=self.silent)):
            self.optimizer.zero_grad()
            loss, mean, var, sigma = self.step(loss_func, n_samples=n_samples, with_loss=True)
            
            accept_epoch = n_max_fail is None or (self.last_n_worse <= n_max_fail)
            
            if accept_epoch:
                if ~(torch.isnan(loss) | torch.isinf(loss)):
                    loss.backward()
                    
                    # Clip exploding gradients
                    torch.nn.utils.clip_grad_norm_(self.nfm.parameters(), 10.)
                    self.optimizer.step()
                    
                    self.call_scheduler(self)
                elif not self.silent:
                    print("Loss is NaN!")
            elif n_max_fail is not None:
                print("StdDev did not improve. Resetting")
                if self.best_params is not None:
                    self.optimizer = torch.optim.Adam(self.nfm.parameters(), **self.optim_args)
                    self.nfm.load_state_dict(self.best_params)
                    self.n_resets += 1
                
            if isinstance(test_callback, Callable):
                test_callback(self)
                
            self.epoch_abs += 1
            # Save summary results            
            res_tot, sigma_tot = variance_weighted_result(np.array(self.statistics['train_results']),
                                                          np.array(self.statistics['train_sigmas']))

            #if i % 10 == 0:
            pbar.set_description(f"Epoch {epoch}: {loss.item():.3f} | σ={sigma:.3f} | σtot={sigma_tot:.3f} | TRes={res_tot:.3f} | MEM free: {self.get_free_mem():.2f} MB")
            
        self.empty_cache()
    
    def evaluate(self, use_adaption:bool=False, use_best:bool=True, n_samples:Optional[int]=None, revert:bool=True):
        if n_samples is not None:
            if use_best and not self.using_best:
                print(f'Loading best params form epoch {self.best_epoch}')
                self.nfm.load_state_dict(self.best_params)
                self.using_best = True
                
                if revert:
                    cparams = self.nfm.state_dict()
            
            with torch.no_grad():
                mean, var, sigma = self.step(n_samples=n_samples, with_loss=False, update_best=False)
                
            if use_best and not self.using_best and revert:
                self.nfm.load_state_dict(cparams)
                
            if not use_adaption:
                return mean.item(), sigma.item()

        if use_adaption:
            means, sigmas = self.statistics['train_results'], self.statistics['train_sigmas']

            return variance_weighted_result(np.array(means), np.array(sigmas))
    
    def sample(self, n_samples:int):
        return self.nfm.sample(n_samples)
            
    def print_model_summary(self, with_structure:bool=False):
        print("Using device " + str(self.device))
        print(f"Model with {sum(p.numel() for p in self.nfm.parameters())} parameters")
        
        if with_structure:
            from torchsummary import summary
            summary(self.nfm)
    
    def save(self, filepath:str):
        torch.save(self.nfm.state_dict(), filepath)
    
    def load(self, filepath:str):
        self.nfm.load_state_dict(torch.load(filepath))
    
    def empty_cache(self):
        torch.cuda.empty_cache()
        gc.collect()
        
    def get_free_mem(self):
        if self.device.type == "cpu":
            return 0
        else:
            return torch.cuda.mem_get_info()[0]/1024/1024 #[1] total memory


# From i-flow
def variance_weighted_result(means:np.ndarray, stddevs:np.ndarray):
    """ Computes weighted mean and stddev of given means and
        stddevs arrays, using Inverse-variance weighting
    """
    assert np.size(means) == np.size(stddevs)
    assert means.shape == stddevs.shape
    variance = 1./np.sum(1./stddevs**2, axis=-1)
    mean = np.sum(means/(stddevs**2), axis=-1)
    mean *= variance
    return mean, np.sqrt(variance)
