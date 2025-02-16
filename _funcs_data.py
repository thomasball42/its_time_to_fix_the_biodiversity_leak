# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 10:53:22 2024

@author: Thomas Ball
"""

import numpy as np

colours_stim = { 
                'Ruminant meat' : "#C90D75",
                'Pig meat'       : "#D64A98",
                'Poultry meat'   : "#D880B1",
                'Dairy'          : "#F7BDDD",
                'Eggs'           : "#FFEDF7",
                
                'Grains'             : "#D55E00",
                "Rice"               : "#D88E53",
                "Soybeans"           : "#DCBA9E",
                
                'Roots and tubers'   : "#0072B2",
                'Vegetables'         : "#4F98C1",
                'Legumes and pulses' : "#9EBFD2",
                
                'Bananas'           : "#FFED00",
                'Tropical fruit'    : "#FFF357",
                'Temperate fruit'   : "#FDF8B9",
                'Tropical nuts'     : "#27E2FF",
                'Temperate nuts'    : "#7DEEFF",
                
                'Sugar beet'    : "#FFC000",
                'Sugar cane'    : "#F7C93B",
                'Spices'        : "#009E73",
                'Coffee'        : "#33CCA2",
                'Cocoa'         : "#62DEBC",
                "Tea and maté"  : "#A2F5DE",
                
                "Oilcrops" : "#000000",
                "Other" : "#A2A2A2"
                }

def weighted_quantile(values, quantiles, sample_weight=None, 
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)