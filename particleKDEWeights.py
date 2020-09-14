#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 17:34:42 2018
Functions to calculate the plastic abundance/mass inside the virtual particles.
@author: kaandorp
"""

import numpy as np

def weight_beachingProbability(t_coastal, pdf_tau):
    """
    calculate probability particle gets removed given a certain time spent near the coast
    pdf_tau: beaching time scale, [days]) 
    """

    weight = np.exp(-t_coastal/pdf_tau)
    weight[weight<0] = 0
    weight[t_coastal==0] = 1

    return weight


def weight_particleSource(source,ratio_popToRiver,ratio_fishToRiver,n_pop,n_river,n_fish,use_kgs=False,kg_particle_river_mid=1.,scale_factor_particle_kg=1.):
    """
    calculate particle weights due to source, i.e. are they coming from rivers or population density
    """
    
    weight = np.zeros(source.shape)

    if use_kgs: #scale directly to kilograms
        
        weight[source == 2] = kg_particle_river_mid
        weight[source == 1] = ratio_popToRiver*kg_particle_river_mid*scale_factor_particle_kg*(n_river/n_pop) #bring to kg's, and scale with the total amount of particles in case one simulation has less

        if n_fish > 0:
            weight[source == 3] = ratio_fishToRiver*kg_particle_river_mid*scale_factor_particle_kg*(n_river/n_fish)
        else:
            weight[source == 3] = 0
            
    else: #scale to river units

        weight[source == 1] = ratio_popToRiver*1.*scale_factor_particle_kg*(n_river/n_pop)
        weight[source == 2] = 1.

        if n_fish > 0:
            weight[source == 3] = ratio_fishToRiver*1.*scale_factor_particle_kg*(n_river/n_fish)
        else:
            weight[source == 3] = 0
 
    return weight


def weight_sinkingBiofouling(t, cdf_tau, cdf_rate, init_sink):
    """
    calculate probability particle gets removed given a certain time age (sinking, biofouling)
    cdf_tau: time scale in weeks
    cdf_rate: sinking rate
    init_sink: fraction of plastics immediately sinking down
    """
    
    rate_min = 1.
    if cdf_rate < rate_min:
        cdf_rate = rate_min
        
    weight = (1-init_sink)*(1 / ( 1 + np.exp((cdf_rate/cdf_tau)*(t-cdf_tau))))
    
    
    weight[weight<0] = 0.
    weight[t==0] = 1.

    return weight    


def weight_river_lowhigh(source,factor,init_month,init_lon_id,init_lat_id,riverInput_low,riverInput_mid,riverInput_high):
    """
    interpolate between low,mid, and high river input estimates from Lebreton et al. (2017)
    """
    if factor > 1:
        factor = 1
    if factor < -1:
        factor = -1
    
    weight = np.ones(source.shape)
    
    riverInput_ratio = np.zeros(riverInput_mid.shape)
    
    indices_sample = np.where(riverInput_mid > 0)
    
    if factor >= 0:
        riverInput_i_v = (1-factor)*riverInput_mid[indices_sample] + factor*riverInput_high[indices_sample]
        riverInput_ratio[indices_sample] = riverInput_i_v / riverInput_mid[indices_sample]
    else:
        riverInput_i_v = (1+factor)*riverInput_mid[indices_sample] - factor*riverInput_low[indices_sample]
        riverInput_ratio[indices_sample] = riverInput_i_v / riverInput_mid[indices_sample]
        
    scale_factor = ((riverInput_ratio*riverInput_mid).sum()) / (riverInput_mid.sum())

    weight[source==2] = riverInput_ratio[init_month[source==2],init_lat_id[source==2],init_lon_id[source==2]]
    
    return weight,scale_factor


def calculateKDEWeights(t_age,t_coastal,source,t_beach_tau,t_sink_tau,t_sink_rate,init_sink,
                        ratio_popToRiver,ratio_fishToRiver,factor_river,init_month,init_lon_id,init_lat_id,riverInput_low,riverInput_mid,riverInput_high,
                        n_pop=1,n_river=1,n_fish=1,kg_particle_river_mid=1.,use_kgs=False):
    """ 
    Calculate particle weights due to all effects taken in account
    """

    t_age_weeks = t_age/ (60*60*24*7)
    t_coastal_days = t_coastal / (60*60*24)


    weight_beaching = weight_beachingProbability(t_coastal_days, t_beach_tau)

    weight_sinking = weight_sinkingBiofouling(t_age_weeks, t_sink_tau, t_sink_rate, init_sink)

    weight_river,scale_factor_particle_kg = weight_river_lowhigh(source,factor_river,init_month,init_lon_id,init_lat_id,riverInput_low,riverInput_mid,riverInput_high)

    weight_source = weight_particleSource(source,ratio_popToRiver,ratio_fishToRiver,n_pop,n_river,n_fish,
                                          use_kgs=use_kgs,kg_particle_river_mid=kg_particle_river_mid,scale_factor_particle_kg=scale_factor_particle_kg)

    weight = weight_beaching * weight_source * weight_sinking * weight_river

    
    return weight

