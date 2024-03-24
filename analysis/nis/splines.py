import torch
from torch.nn import functional as F

import numpy as np

DEFAULT_MIN_BIN_WIDTH = 1e-15
DEFAULT_MIN_BIN_HEIGHT = 1e-15
DEFAULT_MIN_DERIVATIVE = 1e-15

if True:
    def _padded(t:torch.Tensor, lhs):
        lhs = torch.as_tensor(lhs, dtype=t.dtype)
        #zeros = torch.zeros([t.ndim - 1, 2], dtype=torch.int32)
        #lhs_paddings = torch.cat([zeros, torch.tensor([[1, 0]])], axis=0)
        #print('lhs_paddings', lhs_paddings)
        #print('lhs', lhs)
        #print('t', t)
        result = torch.nn.functional.pad(t, pad=(1,0,0,0,0,0), value=lhs)
        #print(result)

        return result

    def _knot_positions(bin_sizes, range_min):
        return _padded(torch.cumsum(bin_sizes, dim=-1, dtype=bin_sizes.dtype) + range_min, lhs=range_min)

    def _gather_squeeze(params, indices):
        rank = indices.dim()
        if rank is None:
            raise ValueError('`indices` must have a statically known rank.')
        
        res = torch.gather(params, dim=-1, index=indices)[..., 0]
        return res

    def _search_sorted(cdf, inputs):
        return torch.maximum(torch.zeros([], dtype=torch.int32), torch.searchsorted(
            cdf[..., :-1], inputs.unsqueeze(-1), side='right'
        ) - 1)
    
    def rational_quadratic_spline(inputs:torch.Tensor, unnormalized_widths:torch.Tensor, unnormalized_heights:torch.Tensor, unnormalized_derivatives:torch.Tensor,
                                inverse:bool=False, left:float=0., right:float=1., bottom:float=0., top:float=1.,
                                min_bin_width:float=DEFAULT_MIN_BIN_WIDTH, min_bin_height:float=DEFAULT_MIN_BIN_HEIGHT,
                                min_derivative:float=DEFAULT_MIN_DERIVATIVE):
        
        #inputs = torch.as_tensor(inputs, dtype=torch.float32)
        #unnormalized_widths = torch.as_tensor(unnormalized_widths, dtype=torch.float32)
        #unnormalized_heights = torch.as_tensor(unnormalized_heights, dtype=torch.float32)
        #unnormalized_derivatives = torch.as_tensor(unnormalized_derivatives, dtype=torch.float32)
        
        out_of_bounds = (inputs < left) | (inputs > right)
        torch.where(out_of_bounds, torch.tensor(left, dtype=inputs.dtype), inputs)

        num_bins = unnormalized_widths.shape[-1]
        # check that number of widths, heights, and derivatives match
        assert num_bins == unnormalized_heights.shape[-1] == unnormalized_derivatives.shape[-1] - 1

        if min_bin_width * num_bins > 1.0:
            raise ValueError('Minimal bin width too large for the number of bins')
        if min_bin_height * num_bins > 1.0:
            raise ValueError('Minimal bin height too large for the number of bins')

        widths = torch.nn.functional.softmax(unnormalized_widths, dim=-1, dtype=inputs.dtype)
        widths = min_bin_width + (1. - min_bin_width * num_bins) * widths
        cumwidths = _knot_positions(widths, 0)
        cumwidths = (right - left) * cumwidths + left
        widths = cumwidths[..., 1:] - cumwidths[..., :-1]

        derivatives = ((min_derivative + torch.nn.functional.softplus(unnormalized_derivatives))
                    / (min_derivative + torch.log(torch.tensor(2., dtype=inputs.dtype))))

        heights = torch.nn.functional.softmax(unnormalized_heights, dim=-1)
        heights = min_bin_height + (1. - min_bin_height * num_bins) * heights
        cumheights = _knot_positions(heights, 0)
        cumheights = (top - bottom) * cumheights + bottom
        heights = cumheights[..., 1:] - cumheights[..., :-1]

        if inverse:
            bin_idx = _search_sorted(cumheights, inputs)
        else:
            bin_idx = _search_sorted(cumwidths, inputs)
        
        #print('bin_idx', bin_idx)
        
        input_cumwidths = _gather_squeeze(cumwidths, bin_idx)
        #print('input_cumwidths', input_cumwidths)
        
        input_bin_widths = _gather_squeeze(widths, bin_idx)

        input_cumheights = _gather_squeeze(cumheights, bin_idx)
        #print('inputs', inputs)
        #print('input_cumheights', input_cumheights)

        delta = heights / widths
        input_delta = _gather_squeeze(delta, bin_idx)
        #print('delta', delta)
        #print('input_delta', input_delta)

        input_derivatives = _gather_squeeze(derivatives, bin_idx)
        input_derivatives_p1 = _gather_squeeze(derivatives[..., 1:], bin_idx)
        #print('input_derivatives', input_derivatives)
        #print('input_derivatives_p1', input_derivatives_p1)

        input_heights = _gather_squeeze(heights, bin_idx)
        #print('input_heights', input_heights)

        if inverse:
            #print('inputs inv', inputs.shape, inputs.dtype, inputs)
            #print('input_cumheights inv', input_cumheights.shape, input_cumheights.dtype, input_cumheights)
            
            aa1 = (inputs - input_cumheights)
            aa2a = input_derivatives + input_derivatives_p1
            aa2b = 2. * input_delta
            aa2 = aa2a - aa2b
            #aa2 = (input_derivatives + input_derivatives_p1 - 2 * input_delta)
            aa = aa1 * aa2
            
            a = aa + input_heights * (input_delta - input_derivatives)
            
            #print('aa2a inv', aa2a.shape, aa2a.dtype, aa2a)
            #print('aa2b inv', aa2b.shape, aa2b.dtype, aa2b)
            #print('aa1 inv', aa1.shape, aa1.dtype, aa1)
            #print('aa2 inv', aa2.shape, aa2.dtype, aa2)
            #print('aa inv', aa.shape, aa.dtype, aa)
            #print('a inv', a.shape, a.dtype, a)
            
            b = (input_heights * input_derivatives
                - (inputs - input_cumheights) * (input_derivatives
                                                + input_derivatives_p1
                                                - 2. * input_delta))
            c = - input_delta * (inputs - input_cumheights)

            discriminant = b ** 2. - 4. * a * c

            theta = (2. * c) / (-b - discriminant.sqrt())
            outputs = theta * input_bin_widths + input_cumwidths

            theta_one_minus_theta = theta * (1. - theta)
            denominator = input_delta + ((input_derivatives + input_derivatives_p1
                                        - 2. * input_delta)
                                        * theta_one_minus_theta)
            derivative_numerator = input_delta ** 2. * (input_derivatives_p1
                                                        * theta ** 2.
                                                        + 2. * input_delta
                                                        * theta_one_minus_theta
                                                        + input_derivatives
                                                        * (1. - theta) ** 2.)
            logabsdet = torch.log(derivative_numerator) - 2. * torch.log(denominator)
            
            #print('logabsdet INV', logabsdet)

            return outputs, -logabsdet
        else:
            theta = (inputs - input_cumwidths) / input_bin_widths
            theta_one_minus_theta = theta * (1. - theta)
            
            #print('theta for', theta)

            numerator = input_heights * (input_delta * theta ** 2.
                                        + input_derivatives
                                        * theta_one_minus_theta)
            denominator = input_delta + ((input_derivatives + input_derivatives_p1
                                        - 2. * input_delta)
                                        * theta_one_minus_theta)
            outputs = input_cumheights + numerator / denominator

            derivative_numerator = input_delta ** 2. * (input_derivatives_p1
                                                        * theta ** 2.
                                                        + 2. * input_delta
                                                        * theta_one_minus_theta
                                                        + input_derivatives
                                                        * (1. - theta) ** 2.)
            logabsdet = torch.log(derivative_numerator) - 2. * torch.log(denominator)
            
            #print('logabsdet FOR', logabsdet)

            return outputs, logabsdet
else:
    def searchsorted(bin_locations, inputs, eps=1e-6):
        bin_locations[..., -1] += eps
        return torch.sum(
            inputs[..., None] >= bin_locations,
            dim=-1
        ) - 1

    def rational_quadratic_spline(inputs,
                                unnormalized_widths,
                                unnormalized_heights,
                                unnormalized_derivatives,
                                inverse=False,
                                left=0., right=1., bottom=0., top=1.,
                                min_bin_width=DEFAULT_MIN_BIN_WIDTH,
                                min_bin_height=DEFAULT_MIN_BIN_HEIGHT,
                                min_derivative=DEFAULT_MIN_DERIVATIVE):
        if torch.min(inputs) < left or torch.max(inputs) > right:
            raise Exception()

        num_bins = unnormalized_widths.shape[-1]

        if min_bin_width * num_bins > 1.0:
            raise ValueError('Minimal bin width too large for the number of bins')
        if min_bin_height * num_bins > 1.0:
            raise ValueError('Minimal bin height too large for the number of bins')

        widths = F.softmax(unnormalized_widths, dim=-1)
        widths = min_bin_width + (1 - min_bin_width * num_bins) * widths
        cumwidths = torch.cumsum(widths, dim=-1)
        cumwidths = F.pad(cumwidths, pad=(1, 0), mode='constant', value=0.0)
        cumwidths = (right - left) * cumwidths + left
        cumwidths[..., 0] = left
        cumwidths[..., -1] = right
        widths = cumwidths[..., 1:] - cumwidths[..., :-1]

        derivatives = min_derivative + F.softplus(unnormalized_derivatives)

        heights = F.softmax(unnormalized_heights, dim=-1)
        heights = min_bin_height + (1 - min_bin_height * num_bins) * heights
        cumheights = torch.cumsum(heights, dim=-1)
        cumheights = F.pad(cumheights, pad=(1, 0), mode='constant', value=0.0)
        cumheights = (top - bottom) * cumheights + bottom
        cumheights[..., 0] = bottom
        cumheights[..., -1] = top
        heights = cumheights[..., 1:] - cumheights[..., :-1]

        if inverse:
            bin_idx = searchsorted(cumheights, inputs)[..., None]
        else:
            bin_idx = searchsorted(cumwidths, inputs)[..., None]

        input_cumwidths = cumwidths.gather(-1, bin_idx)[..., 0]
        input_bin_widths = widths.gather(-1, bin_idx)[..., 0]

        input_cumheights = cumheights.gather(-1, bin_idx)[..., 0]
        delta = heights / widths
        input_delta = delta.gather(-1, bin_idx)[..., 0]

        input_derivatives = derivatives.gather(-1, bin_idx)[..., 0]
        input_derivatives_plus_one = derivatives[..., 1:].gather(-1, bin_idx)[..., 0]

        input_heights = heights.gather(-1, bin_idx)[..., 0]

        if inverse:
            a = (((inputs - input_cumheights) * (input_derivatives
                                                + input_derivatives_plus_one
                                                - 2 * input_delta)
                + input_heights * (input_delta - input_derivatives)))
            b = (input_heights * input_derivatives
                - (inputs - input_cumheights) * (input_derivatives
                                                + input_derivatives_plus_one
                                                - 2 * input_delta))
            c = - input_delta * (inputs - input_cumheights)

            discriminant = b.pow(2) - 4 * a * c
            assert (discriminant >= 0).all()

            root = (2 * c) / (-b - torch.sqrt(discriminant))
            outputs = root * input_bin_widths + input_cumwidths

            theta_one_minus_theta = root * (1 - root)
            denominator = input_delta + ((input_derivatives + input_derivatives_plus_one - 2 * input_delta)
                                        * theta_one_minus_theta)
            derivative_numerator = input_delta.pow(2) * (input_derivatives_plus_one * root.pow(2)
                                                        + 2 * input_delta * theta_one_minus_theta
                                                        + input_derivatives * (1 - root).pow(2))
            logabsdet = torch.log(derivative_numerator) - 2 * torch.log(denominator)

            return outputs, -logabsdet
        else:
            theta = (inputs - input_cumwidths) / input_bin_widths
            theta_one_minus_theta = theta * (1 - theta)

            numerator = input_heights * (input_delta * theta.pow(2)
                                        + input_derivatives * theta_one_minus_theta)
            denominator = input_delta + ((input_derivatives + input_derivatives_plus_one - 2 * input_delta)
                                        * theta_one_minus_theta)
            outputs = input_cumheights + numerator / denominator

            derivative_numerator = input_delta.pow(2) * (input_derivatives_plus_one * theta.pow(2)
                                                        + 2 * input_delta * theta_one_minus_theta
                                                        + input_derivatives * (1 - theta).pow(2))
            logabsdet = torch.log(derivative_numerator) - 2 * torch.log(denominator)

            return outputs, logabsdet