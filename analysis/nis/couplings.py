import torch
import torch.nn as nn
import numpy as np
from normflows.flows import Flow
from analysis.nis.splines import rational_quadratic_spline, DEFAULT_MIN_BIN_HEIGHT, DEFAULT_MIN_BIN_WIDTH, DEFAULT_MIN_DERIVATIVE
from typing import Union

class CouplingBijector(Flow):
    def __init__(self, mask, transform_net_create_fn, options=None, **kwargs):
        mask = torch.as_tensor(mask)
        super().__init__()

        self.features = len(mask)
        features_vector = torch.arange(self.features)
        
        self.register_buffer('identity_features', tensor=features_vector.masked_select(mask <= 0))
        self.register_buffer('transform_features', tensor=features_vector.masked_select(mask > 0))
        
        #print('Init with identity_features', self.identity_features, 'and transform_features', self.transform_features)
        
        assert (self.num_identity_features + self.num_transform_features == self.features)

        self.transform_net = transform_net_create_fn(
            self.num_identity_features,
            self.num_transform_features * self._transform_dim_multiplier(),
            options
        )

    @property
    def num_identity_features(self):
        # Number of un-transformed features
        return len(self.identity_features)

    @property
    def num_transform_features(self):
        # Number of transformed features
        return len(self.transform_features)

    def forward(self, inputs, context=None):
        identity_split = inputs[:, self.identity_features]
        transform_split = inputs[:, self.transform_features]

        transform_params = self.transform_net(identity_split, context)
            
        transform_split, logabsdet = self._coupling_transform_forward(
            inputs=transform_split,
            transform_params=transform_params
        )

        outputs = torch.cat([identity_split, transform_split], dim=1)
        indices = torch.cat([self.identity_features, self.transform_features], dim=-1)
        outputs = outputs[:, torch.argsort(indices)]
        
        #print("indices", torch.argsort(indices))
        
        return outputs, logabsdet

    def inverse(self, inputs, context=None):
        identity_split = inputs[:, self.identity_features]
        transform_split = inputs[:, self.transform_features]

        transform_params = self.transform_net(identity_split, context)

        transform_split, logabsdet = self._coupling_transform_inverse(
            inputs=transform_split,
            transform_params=transform_params
        )

        outputs = torch.cat([identity_split, transform_split], dim=1)
        indices = torch.cat([self.identity_features, self.transform_features], dim=-1)
        outputs = outputs[:, torch.argsort(indices)]

        return outputs, logabsdet

    def _forward_log_det_jacobian(self, inputs, context=None):
        """ Compute forward log det Jacobian. """
        identity_split = inputs[:, self.identity_features]
        transform_split = inputs[:, self.transform_features]

        transform_params = self.transform_net(identity_split, context)

        transform_split, logabsdet = self._coupling_transform_forward(
            inputs=transform_split,
            transform_params=transform_params
        )

        return logabsdet

    def _inverse_log_det_jacobian(self, inputs, context=None):
        """ Compute Inverse log det Jacobian. """

        identity_split = inputs[:, self.identity_features]
        transform_split = inputs[:, self.transform_features]

        transform_params = self.transform_net(identity_split, context)

        transform_split, logabsdet = self._coupling_transform_inverse(
            inputs=transform_split,
            transform_params=transform_params
        )

        return logabsdet

    def _transform_dim_multiplier(self):
        raise NotImplementedError()

    def _coupling_transform_forward(self, inputs, transform_params):
        raise NotImplementedError()

    def _coupling_transform_inverse(self, inputs, transform_params):
        raise NotImplementedError()

class PiecewiseBijector(CouplingBijector):
    def _coupling_transform_forward(self, inputs:torch.Tensor, transform_params:torch.Tensor):
        #print(f'PiecewiseBijector[{self.i_test}] : FOR')
        return self._coupling_transform(inputs, transform_params, inverse=False)

    def _coupling_transform_inverse(self, inputs:torch.Tensor, transform_params:torch.Tensor):
        #print(f'PiecewiseBijector[{self.i_test}] : INV')
        return self._coupling_transform(inputs, transform_params, inverse=True)

    def _coupling_transform(self, inputs:torch.Tensor, transform_params:torch.Tensor, inverse=False):
        b, d = inputs.shape
        transform_params = transform_params.view(b, d, -1)

        outputs, logabsdet = self._piecewise_cdf(inputs, transform_params, inverse)
        
        logabsdet = torch.sum(logabsdet, dim=-1)
        
        #print(f'{"INV" if inverse else "FOR"} [{self.i_test}]', logabsdet)

        return outputs, logabsdet

    def _piecewise_cdf(self, inputs, transform_params, inverse=False)->tuple[torch.Tensor,torch.Tensor]:
        raise NotImplementedError()

class PiecewiseRationalQuadraticSpline(PiecewiseBijector):
    def __init__(self, mask, transform_net_create_fn, num_bins=10,
                 min_bin_width=DEFAULT_MIN_BIN_WIDTH,
                 min_bin_height=DEFAULT_MIN_BIN_HEIGHT,
                 min_derivative=DEFAULT_MIN_DERIVATIVE,
                 i_test=0,
                 **kwargs):
        self.num_bins = num_bins
        self.min_bin_width = min_bin_width
        self.min_bin_height = min_bin_height
        self.min_derivative = min_derivative

        super(PiecewiseRationalQuadraticSpline, self).__init__(
            mask, transform_net_create_fn, **kwargs)
        
        self.i_test = i_test

    def _transform_dim_multiplier(self):
        return self.num_bins * 3 + 1

    def _piecewise_cdf(self, inputs, transform_params, inverse=False):
        unnormalized_widths = transform_params[..., :self.num_bins]
        unnormalized_heights = transform_params[...,
                                                self.num_bins:2*self.num_bins]
        unnormalized_derivatives = transform_params[..., 2*self.num_bins:]

        res = rational_quadratic_spline(
            inputs=inputs,
            unnormalized_widths=unnormalized_widths,
            unnormalized_heights=unnormalized_heights,
            unnormalized_derivatives=unnormalized_derivatives,
            inverse=inverse,
            min_bin_width=self.min_bin_width,
            min_bin_height=self.min_bin_height,
            min_derivative=self.min_derivative
        )
        
        #print(f'COUPL {self.i_test}')
        
        return res