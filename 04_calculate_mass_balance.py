#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 12:21:53 2019
After finding the optimal set of parameters, we can now estimate in which environmental compartments we have all the plastic map. Results are plotted.

@author: kaandorp
"""

from datetime import datetime, timedelta
import numpy as np
import glob
from netCDF4 import Dataset
import pandas as pd
import os
import particleKDEWeights
import pickle

def getLandBorder(landMask,lon,lat,val_add):
    """
    Function to obtain a mask of the land which borders ocrean, uses the landmask and searches for 
    boundaries with the sea (horiz./vert. adjacent cells only)
    TODO: check for cyclic boundaries
    TODO: check diagonal cells as well?
    """
    
    n_lat = landMask.shape[0]
    n_lon = landMask.shape[1]

    for i1 in range(n_lat):
        for i2 in range(n_lon):
            
            check_bot = True
            check_top = True
            check_left = True
            check_right = True
            
            # check whether land is located at boundary
            if i1 == 0:
                check_top = False
            if i1 == n_lat-1:
                check_bot = False
            if i2 == 0:
                check_left = False
            if i2 == n_lon-1:
                check_right = False
                
            # check whether cell is land, if so look for coast
            if landMask[i1,i2] == 1:
                
                if check_top:
                    if (landMask[i1-1,i2] == 0) or (landMask[i1-1,i2] >= 2):
                        landMask[i1,i2] = -1
                if check_bot:
                    if (landMask[i1+1,i2] == 0) or (landMask[i1+1,i2] >= 2):
                        landMask[i1,i2] = -1
                if check_left:
                    if (landMask[i1,i2-1] == 0) or (landMask[i1,i2-1] >= 2):
                        landMask[i1,i2] = -1
                if check_right:
                    if (landMask[i1,i2+1] == 0) or (landMask[i1,i2+1] >= 2):
                        landMask[i1,i2] = -1
    landMask[landMask == -1] = val_add
            
    return landMask
#%% 
lons = np.linspace(-6,36.25,677)
lats = np.linspace(30.1875,45.9375,253)

dataDir = 'Data/Temp/34_delete180Days/K10'

paramFile = './result_files/Results_180Days_K10.pickle'

outFile = './result_files/hist_180Days_K10.pickle'


Kval = 10
time_origin = datetime(2006,2,1,0,0)
n_folders_per_source = 6 #limit computational costs if necessary

use_pop = True
use_rivers = True
use_fisheries = True
PLOT = True


use_SSD = True

if os.environ['USER'] == 'kaandorp': # desktop
    
    if use_SSD:
        homeDir = '/Volumes/externe_SSD/kaandorp/'
    else:
        homeDir = os.environ['HOME']
        
    all_folders = []
    
    if use_pop:
        particleFolders_pop = os.path.join(homeDir , dataDir, 'popOnly/')
        folders_pop = np.sort(glob.glob(os.path.join(particleFolders_pop,'*')))[0:n_folders_per_source]
        all_folders.append(folders_pop)
    if use_rivers:
        particleFolders_rivers = os.path.join(homeDir , dataDir, 'riversOnly/')
        folders_rivers = np.sort(glob.glob(os.path.join(particleFolders_rivers,'*')))[0:n_folders_per_source]   
        all_folders.append(folders_rivers)
    if use_fisheries:
        particleFolders_fisheries = os.path.join(homeDir , dataDir, 'fisheriesOnly/')
        folders_fisheries = np.sort(glob.glob(os.path.join(particleFolders_fisheries,'*')))[0:n_folders_per_source]   
        all_folders.append(folders_fisheries)
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import cmocean


    
 
#%%
lons_grid = np.linspace(-6,36.25,677)
lats_grid = np.linspace(30.1875,45.9375,253)
fieldMesh_x,fieldMesh_y = np.meshgrid(lons_grid,lats_grid)

lon_spacing = lons_grid[1] - lons_grid[0]
lat_spacing = lats_grid[1] - lats_grid[0]
lons_grid2 = lons_grid - 0.5*lon_spacing
lats_grid2 = lats_grid - 0.5*lat_spacing
lons_grid2 = np.append(lons_grid2,lons_grid2[-1]+lon_spacing)
lats_grid2 = np.append(lats_grid2,lats_grid2[-1]+lat_spacing)

fieldMesh_xA, fieldMesh_yA = np.meshgrid(lons_grid2,lats_grid2)

if not os.path.exists(outFile):  
    folder_0_pop = folders_pop[0]
    fileList_p0_pop = list(np.sort(glob.glob(os.path.join(folder_0_pop ,'K%3.1f_data*' % Kval))))
        
    files_initStates = [glob.glob(os.path.join(folder,'particlesRelease*')) for folder in np.concatenate((folders_pop,folders_rivers,folders_fisheries))]
    list_ = []
    
    
    index_plus = [0] #since particle files are concatenated, indices between different files need to kept track of
    for i,file in enumerate(files_initStates):
        df_i = pd.read_csv(file[0],header=None,sep='\t')
        df_i.set_index(df_i.index+index_plus[i], inplace=True)
        list_.append(df_i)
        
        index_plus.append(index_plus[i]+len(df_i))
    
    index_plus= np.array(index_plus)    
    initStates = pd.concat(list_)
    init_indices = initStates.index
    init_time = initStates.iloc[:,0]
    init_lat = initStates.iloc[:,1]
    init_lon = initStates.iloc[:,2]
    init_popDen = initStates.iloc[:,3]
    init_source = np.zeros(init_popDen.shape)
    init_source[0:index_plus[len(folders_pop)]] = 1
    init_source[index_plus[len(folders_pop)]:index_plus[len(folders_pop)+len(folders_rivers)]] = 2
    init_source[index_plus[len(folders_pop)+len(folders_rivers)]: ] = 3
    
    init_lon_id = np.digitize(init_lon,lons_grid,right=True)
    init_lat_id = np.digitize(init_lat,lats_grid,right=True)
    
    
    time_20100101 = (datetime(2010,1,1,0,0)-time_origin).total_seconds()
    time_20101231 = (datetime(2010,12,31,23,59)-time_origin).total_seconds()
    
    n_pop = len(init_source[init_source == 1])
    n_river = len(init_source[init_source == 2])
    n_fish = len(init_source[init_source == 3])
    
    n_river_2010 = len(init_source[(init_source == 2) & (init_time >= time_20100101) & (init_time <= time_20101231)])
    kg_year_rivers_mid = 1248000 # from Lebreton (2017), mid estimate
    kg_particle_river_mid = kg_year_rivers_mid / n_river_2010 
    
    init_month_id = np.array([ (time_origin + timedelta(seconds=int(init_time[i1]))).month - 1 for i1 in range(len(init_time.values))])
    riverInput_low = np.loadtxt('datafile_riverInputMatrix_677x_253y_low').reshape([12,len(lats),len(lons)])
    riverInput_mid = np.loadtxt('datafile_riverInputMatrix_677x_253y').reshape([12,len(lats),len(lons)])
    riverInput_high = np.loadtxt('datafile_riverInputMatrix_677x_253y_high').reshape([12,len(lats),len(lons)])
    
    vals_labels = ['beach_tau','sink_tau','sink_rate','sink_init','ratio_popToRiver','ratio_tourToRiver','ratio_fishToRiver','KDE_bw','river_imp'] 
    
    
    pickle_file = open(paramFile,'rb')
    dict_results = pickle.load(pickle_file)
        
    def to_unitfull(x_unitless_,x_lower_,x_upper_,log_scale_):
        vals_unitfull_ = np.zeros(x_lower_.shape)
        for i1,bool_log in enumerate(log_scale_):
            if bool_log:
                vals_unitfull_[i1] = 10**(x_lower_[i1]+x_unitless_[i1]*(x_upper_[i1]-x_lower_[i1]))
            else:
                vals_unitfull_[i1] = x_lower_[i1]+x_unitless_[i1]*(x_upper_[i1]-x_lower_[i1])
        return vals_unitfull_
    
    x_op = to_unitfull(dict_results['x_op_dimless'],dict_results['x_lower'],dict_results['x_upper'],dict_results['log_scale'])
    
    
    t_beach_tau = x_op[0]
    t_sink_tau = x_op[1]
    t_sink_rate = x_op[2]
    sink_init = x_op[3]
    ratio_popToRiver = x_op[4]
    ratio_fishToRiver = x_op[5]
    KDE_bw = x_op[6]
    factor_river = x_op[7]
    
      
    type_deg = 'weight'
    use_kgs = True  
    
    ##
    n_months = 11*12-1
    mass_floating = np.zeros([14*len(fileList_p0_pop)])
    mass_cons_acc = np.zeros([14*len(fileList_p0_pop)]) #mass conservation accuracy
    
    i0 = 0 # index counting : i3 + 3*i2 ?
    list_masses_beach = np.zeros([index_plus[-1],14*len(fileList_p0_pop)])
    list_masses_sink = np.zeros([index_plus[-1],14*len(fileList_p0_pop)])
    
    index_list_masses_col = 0
    
    old_date = datetime(2006,2,1)
    weight_added_pop = np.zeros([14*len(fileList_p0_pop)])
    weight_added_riv = np.zeros([14*len(fileList_p0_pop)])
    weight_added_fis = np.zeros([14*len(fileList_p0_pop)])
    weight_lost_sink = np.zeros([14*len(fileList_p0_pop)])
    weight_lost_beach = np.zeros([14*len(fileList_p0_pop)])
    
    print('Calculating histograms...')
    
    hist_beaching = np.zeros([n_months,len(lons),len(lats)])
    hist_sinking = np.zeros([n_months,len(lons),len(lats)])
    hist_floating = np.zeros([n_months,len(lons),len(lats)])
    weight_tot = np.zeros([14*len(fileList_p0_pop)])
    traj_seen = np.array([])
    
    for i1 in range(len(fileList_p0_pop)):
    
        counter_i2i3 = 0
    
        # go through all folders per type
        for i2,use_type in enumerate([use_pop,use_rivers,use_fisheries]):
            
            if use_type:
                
                for i3 in range(len(all_folders[i2])): #3 folders per type
                    tmp_filelist = list(np.sort(glob.glob(os.path.join(all_folders[i2][i3] ,'K%3.1f_data*' % Kval))))
                    tmp_data = Dataset(tmp_filelist[i1])
                    
                                    

                    if (i2 == 0) and (i3 == 0):
                    
                        t_age = tmp_data['age'][:]
                        t_coastal = tmp_data['coastalAge_1cell'][:]
                        lon = tmp_data['lon'][:]
                        lat = tmp_data['lat'][:]
                        traj_ = tmp_data['trajectory'][:,0] + index_plus[counter_i2i3]
                        source = (i2+1)*np.ones(traj_.shape)
                        init_month_id_ = init_month_id[traj_]
                        init_lon_id_ = init_lon_id[traj_]
                        init_lat_id_ = init_lat_id[traj_]
                        time_ = tmp_data['time'][:]
                        init_index = traj_.copy()
                        
                        len_time = len(time_[:,0])
        
                    else:
                        t_age = np.vstack((t_age,tmp_data['age'][:]))
                        lon = np.vstack((lon,tmp_data['lon'][:]))
                        lat = np.vstack((lat,tmp_data['lat'][:]))                    
                        t_coastal = np.vstack((t_coastal,tmp_data['coastalAge_1cell'][:]))
                        traj_ = tmp_data['trajectory'][:,0] + index_plus[counter_i2i3]
                        source = np.append(source,(i2+1)*np.ones(traj_.shape))
                        init_month_id_ = np.append(init_month_id_,init_month_id[traj_])
                        init_lon_id_ = np.append(init_lon_id_,init_lon_id[traj_])
                        init_lat_id_ = np.append(init_lat_id_,init_lat_id[traj_])
                        time_ = np.vstack((time_,tmp_data['time'][:]))
                        
                        init_index = np.append(init_index,traj_)
                    
                    counter_i2i3 += 1
    
    
        for i4 in range(1,time_.shape[1]):                       
    
            if i1 < 5:
                len_samples = int(1.2*(len_time/15))
                values, counts = np.unique(time_[0:len_samples,i4][~time_[0:len_samples,i4].mask], return_counts=True) 
                curr_time = values[np.argmax(counts)]
            elif i1 < 30:
                values, counts = np.unique(time_[:,i4][~time_[:,i4].mask], return_counts=True) 
                curr_time = values[np.argmax(counts)]
            else:
                curr_time = np.nanmedian(time_[:,i4])  #save computational costs            
                
            
            
            curr_date = datetime(2001,1,1) + timedelta(seconds = curr_time)
            year_ = curr_date.year
            month_ = curr_date.month
            
            if (curr_date - old_date).days != 1:
                print('------------WARNING--------------------')
                print('Days in analysis are not consecutive')
    
                
            selection_ = np.where(time_ == curr_time)
                    
            selection_row = selection_[0][selection_[1] != 0] #remove selection entries which only have the last time (in this case we cannot compute difference)
            selection_col = selection_[1][selection_[1] != 0]
            
            index_list_masses_row = init_index[selection_row] # to check mass conservation, all masses are stored in list_masses. These are the corresponding rows
            
            
            t_age_days = t_age / (60*60*24)
            t_age_weeks = t_age/ (60*60*24*7)
            t_coastal_days = t_coastal / (60*60*24)
           
            weight_river,scale_factor_particle_kg = particleKDEWeights.weight_river_lowhigh(source[selection_row],factor_river,init_month_id_[selection_row],init_lon_id_[selection_row],init_lat_id_[selection_row],riverInput_low,riverInput_mid,riverInput_high)
        
            weight_source = particleKDEWeights.weight_particleSource(source[selection_row],ratio_popToRiver,ratio_fishToRiver,n_pop,n_river,n_fish,
                                                  use_kgs=use_kgs,kg_particle_river_mid=kg_particle_river_mid,scale_factor_particle_kg=scale_factor_particle_kg)
        
            weight_at_source = weight_source*weight_river
    
            fraction_beaching_current   = particleKDEWeights.weight_beachingProbability(t_coastal_days[selection_row,selection_col], t_beach_tau)
            fraction_beaching_previous  = particleKDEWeights.weight_beachingProbability(t_coastal_days[selection_row,selection_col-1], t_beach_tau)
        
            fraction_sinking_current    = particleKDEWeights.weight_sinkingBiofouling(t_age_weeks[selection_row,selection_col], t_sink_tau, t_sink_rate, sink_init)
            fraction_sinking_previous   = particleKDEWeights.weight_sinkingBiofouling(t_age_weeks[selection_row,selection_col-1], t_sink_tau, t_sink_rate, sink_init)
            
            d_fraction_beaching = (fraction_beaching_previous-fraction_beaching_current)
            d_fraction_sinking = (fraction_sinking_previous-fraction_sinking_current)
            
            mass_beaching_lost = np.array(weight_at_source*fraction_sinking_previous*d_fraction_beaching - weight_at_source*(d_fraction_sinking / (d_fraction_sinking+d_fraction_beaching))*d_fraction_sinking*d_fraction_beaching)
            mass_sinking_lost = np.array(weight_at_source*fraction_beaching_previous*d_fraction_sinking - weight_at_source*(d_fraction_beaching / (d_fraction_sinking+d_fraction_beaching))*d_fraction_sinking*d_fraction_beaching)
            mass_floating_remaining = weight_at_source*(fraction_beaching_current*fraction_sinking_current)
        
            weight_added_pop[index_list_masses_col] = weight_at_source[(t_age_days[selection_row,selection_col-1] == 0) & (source[selection_row] == 1)].sum()
            weight_added_riv[index_list_masses_col] = weight_at_source[(t_age_days[selection_row,selection_col-1] == 0) & (source[selection_row] == 2)].sum()
            weight_added_fis[index_list_masses_col] = weight_at_source[(t_age_days[selection_row,selection_col-1] == 0) & (source[selection_row] == 3)].sum()
            
            lon_mean = .5*(lon[selection_row,selection_col] + lon[selection_row,selection_col-1])
            lat_mean = .5*(lat[selection_row,selection_col] + lat[selection_row,selection_col-1])
        
            hist_id_lon = np.digitize(lon_mean,lons_grid,right=True)-1
            hist_id_lat = np.digitize(lat_mean,lats_grid,right=True)-1
              
            index_mat = (year_ - 2006)*12 + (month_-2)
      
            for i_lo,i_la,m_b,m_s,m_f in zip(hist_id_lon,hist_id_lat,mass_beaching_lost,mass_sinking_lost,mass_floating_remaining):
                if ~np.isnan(m_b):
                    hist_beaching[index_mat,i_lo,i_la] += m_b
                if ~np.isnan(m_s):
                    hist_sinking[index_mat,i_lo,i_la] += m_s
                if ~np.isnan(m_f):
                    hist_floating[index_mat,i_lo,i_la] += m_f
                    
            weight_lost_sink[index_list_masses_col] = mass_sinking_lost.sum()
            weight_lost_beach[index_list_masses_col] = mass_beaching_lost.sum()
        
            mass_floating[index_list_masses_col] = np.sum(weight_at_source*(fraction_beaching_current*fraction_sinking_current))
        
            mass_cons_acc[index_list_masses_col] = (weight_lost_beach.sum()+weight_lost_sink.sum()+mass_floating[index_list_masses_col]) / np.sum(weight_tot)
         
            index_list_masses_col += 1
            old_date = curr_date
      
        print(i1)
        
    dict_results = {}
    dict_results['x'] = x_op.copy()
    dict_results['hist_beaching'] = hist_beaching.copy()
    dict_results['hist_sinking'] = hist_sinking.copy()
    dict_results['hist_floating'] = hist_floating.copy()
    dict_results['mass_floating'] = mass_floating
    dict_results['weight_added_pop'] = weight_added_pop
    dict_results['weight_added_riv'] = weight_added_riv
    dict_results['weight_added_fis'] = weight_added_fis
    
    outfile = open(outFile,'wb')
    pickle.dump(dict_results,outfile)
    outfile.close()
        
    weight_added_total = weight_added_pop.sum() + weight_added_riv.sum() + weight_added_fis.sum()
    weight_lost_total =   hist_beaching.sum() + hist_sinking.sum() 
    weight_compartments_total = weight_lost_total + mass_floating[index_list_masses_col-1]

else:
    print('loading histograms...')
    with open(outFile, 'rb') as buff:
        data = pickle.load(buff)    
    hist_beaching = data['hist_beaching']
    hist_sinking = data['hist_sinking']
    hist_floating = data['hist_floating']
    mass_floating = data['mass_floating']
    weight_added_pop = data['weight_added_pop']
    weight_added_riv = data['weight_added_riv']
    weight_added_fis  = data['weight_added_fis']

    weight_added_total = weight_added_pop.sum() + weight_added_riv.sum() + weight_added_fis.sum()
    weight_lost_total =   hist_beaching.sum() + hist_sinking.sum() 
    weight_compartments_total = weight_lost_total + mass_floating[-1]

#%% mass overview
    
print(outFile)
print('percentage sunk: %f' % (hist_sinking.sum()/weight_compartments_total ))
print('percentage beached: %f' % (hist_beaching.sum()/weight_compartments_total ))
print('percentage floating: %f' % (mass_floating[-1]/weight_compartments_total ))
print('mass floating: %f' % (mass_floating[-1]))


if PLOT:
    
    date_arr = pd.date_range(datetime(2006,2,1),datetime(2016,12,20))
    fig,ax = plt.subplots(1)
    l1 = ax.plot_date(date_arr,mass_floating*2,'k-',linewidth=1.5,zorder=1,label='floating plastic mass')
    ax2 = plt.twinx(ax)
    l2 = ax2.plot_date(date_arr,weight_added_pop,'-',linewidth=0.5,zorder=2,label='added mass, coastal population')
    l3 = ax2.plot_date(date_arr,weight_added_riv,'-',linewidth=0.5,zorder=2,label='added mass, riverine')
    l4 = ax2.plot_date(date_arr,weight_added_fis,'-',linewidth=0.5,zorder=2,label='added mass, fisheries')
    
    lns = l1+l2+l3+l4
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs,loc='upper left',fontsize=12)
    
    ax2.set_ylim(0,10000)
    ax.set_xlabel('Date',fontsize=14)
    ax.set_ylabel('Floating plastic mass [kg]',fontsize=14)
    ax2.set_ylabel('Added mass [kg/day]',fontsize=14)
    
    
    i_2015 = (date_arr >= datetime(2015,1,1)) & (date_arr <= datetime(2015,12,31))
    
    weight_added_2015 = weight_added_pop[i_2015].sum() + weight_added_riv[i_2015].sum() + weight_added_fis[i_2015].sum()
    
    #%%
    hist_beaching_sum = hist_beaching.sum(axis=0).T
    hist_sinking_sum = hist_sinking.sum(axis=0).T
    hist_floating_sum = hist_floating.sum(axis=0).T
    
    
    grid_area = np.zeros(fieldMesh_y.shape)
    for i1 in range(grid_area.shape[0]):
        grid_area[i1,:] = (lat_spacing*1.11e2) * (lon_spacing*np.cos(lats[i1]*(np.pi/180))*1.11e2)
            
    #%%
    import calendar
    
    PLOT_MONTHS = False
    
    for i1 in range(12):
        month_arr = pd.date_range(datetime(2006,2,1),datetime(2016,12,1),freq='MS')

        
        log_min = -2
        log_max = 2
        
        arr_id_start = np.roll(np.arange(0,12),1)[i1]
        
        month_arr_ind = np.arange(arr_id_start,131,12)
        
        index_2015_month = 107+i1
        
        if PLOT_MONTHS:
            plt.figure(figsize=(20,6.5))     
            ax = plt.axes(projection=ccrs.PlateCarree())
            ax.coastlines(resolution='10m',color='grey')
            cplot = ax.contourf(fieldMesh_x,fieldMesh_y,hist_floating[month_arr_ind,:,:].mean(axis=0).T/grid_area ,cmap=cmocean.cm.matter,levels=np.logspace(log_min,log_max,100),norm=colors.LogNorm(),transform=ccrs.PlateCarree(),extend='both') #,norm=colors.LogNorm(vmin=1, vmax=500),cmap=plt.cm.OrRd,transform=ccrs.PlateCarree())
            cbar = plt.colorbar(cplot,ticks=[10**i for i in np.linspace(log_min,log_max,log_max-log_min+1)])
            cbar.ax.tick_params(labelsize='x-large')
            cbar.set_label('floating, %s over 2006-2016 [kg/km$^2$]' % calendar.month_name[i1+1],fontsize='xx-large')
    
        n_days = calendar.monthrange(month_arr[index_2015_month].year,month_arr[index_2015_month].month)[1]
        hist_floating_month_2015_avg = hist_floating[index_2015_month,:,:] / n_days
    
        mask_nonzero = (hist_floating_month_2015_avg != 0)
    
        log_floating = np.log10(hist_floating_month_2015_avg[mask_nonzero])
        sigma_meas = np.sqrt(0.2201)
                

        hist_perturbed = np.zeros(hist_floating_month_2015_avg.shape)
        n_times = 500
        for i2 in range(n_times):
            log_floating_perturbed = log_floating +  np.random.normal(0,sigma_meas,len(log_floating))
            hist_perturbed_temp =     hist_floating_month_2015_avg.copy()
            hist_perturbed_temp[mask_nonzero] = 10**log_floating_perturbed
            
            hist_perturbed += hist_perturbed_temp
        
        hist_perturbed /= n_times
        print(hist_floating_month_2015_avg.sum(),hist_perturbed.sum())
    
        if PLOT_MONTHS:
        
            plt.figure(figsize=(20,6.5))     
            ax = plt.axes(projection=ccrs.PlateCarree())
            ax.coastlines(resolution='10m',color='grey')
            cplot = ax.contourf(fieldMesh_x,fieldMesh_y,hist_floating_month_2015_avg.T ,cmap=cmocean.cm.matter,levels=np.logspace(log_min,log_max,100),norm=colors.LogNorm(),transform=ccrs.PlateCarree(),extend='both') #,norm=colors.LogNorm(vmin=1, vmax=500),cmap=plt.cm.OrRd,transform=ccrs.PlateCarree())
            cbar = plt.colorbar(cplot,ticks=[10**i for i in np.linspace(log_min,log_max,log_max-log_min+1)])
            cbar.ax.tick_params(labelsize='x-large')
            cbar.set_label('floating, %s over 2006-2016 [kg/km$^2$]' % calendar.month_name[i1+1],fontsize='xx-large')    
            
            plt.figure(figsize=(20,6.5))     
            ax = plt.axes(projection=ccrs.PlateCarree())
            ax.coastlines(resolution='10m',color='grey')
            cplot = ax.contourf(fieldMesh_x,fieldMesh_y,hist_perturbed.T ,cmap=cmocean.cm.matter,levels=np.logspace(log_min,log_max,100),norm=colors.LogNorm(),transform=ccrs.PlateCarree(),extend='both') #,norm=colors.LogNorm(vmin=1, vmax=500),cmap=plt.cm.OrRd,transform=ccrs.PlateCarree())
            cbar = plt.colorbar(cplot,ticks=[10**i for i in np.linspace(log_min,log_max,log_max-log_min+1)])
            cbar.ax.tick_params(labelsize='x-large')
            cbar.set_label('floating, %s over 2006-2016 [kg/km$^2$]' % calendar.month_name[i1+1],fontsize='xx-large')    
                
        
    
    #%% beaching
    
    log_min = -1
    log_max = 2
    fig = plt.figure(figsize=(3.4,2.2),dpi=120)     
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(resolution='10m',color='grey',linewidth=0.5,zorder=0)
    cplot = ax.contourf(fieldMesh_x,fieldMesh_y,hist_beaching_sum/grid_area ,cmap=cmocean.cm.matter,levels=np.logspace(log_min,log_max,100),norm=colors.LogNorm(),transform=ccrs.PlateCarree(),extend='both',zorder=1) #,norm=colors.LogNorm(vmin=1, vmax=500),cmap=plt.cm.OrRd,transform=ccrs.PlateCarree())
    cax = fig.add_axes([0.15, 0.22, 0.7, 0.02])
    cbar = plt.colorbar(cplot,cax=cax,orientation='horizontal',ticks=[10**i for i in np.linspace(log_min,log_max,log_max-log_min+1)])
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label('beaching, 2006-2016 [kg/km$^2$]',fontsize=7)
    
    
    
    #%%
    # sinking

    log_min = -5
    log_max = 2
    fig = plt.figure(figsize=(3.4,2.2),dpi=120)     
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(resolution='10m',color='grey',linewidth=0.5,zorder=0)
    #cplot = ax.contourf(fieldMesh_x,fieldMesh_y,hist_beaching_sum.sum(axis=0).T ,cmap=plt.cm.OrRd,transform=ccrs.PlateCarree(),norm=colors.LogNorm(vmin=.1, vmax=hist_beaching_sum.sum(axis=0).max())) #,norm=colors.LogNorm(vmin=1, vmax=500),cmap=plt.cm.OrRd,transform=ccrs.PlateCarree())
    cplot = ax.contourf(fieldMesh_x,fieldMesh_y,hist_sinking_sum / grid_area,levels=np.logspace(log_min,log_max,100),cmap=cmocean.cm.matter,norm=colors.LogNorm(),transform=ccrs.PlateCarree(),extend='both') #,norm=colors.LogNorm(vmin=1, vmax=500),cmap=plt.cm.OrRd,transform=ccrs.PlateCarree())
    cax = fig.add_axes([0.15, 0.22, 0.7, 0.02])
    
    cbar = plt.colorbar(cplot,cax=cax,orientation='horizontal',ticks=[10**i for i in np.linspace(log_min,log_max,log_max-log_min+1)])
    
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label('sinking, 2006-2016 [kg/km$^2$]',fontsize=7)
    
    
    #%%
    
    #def closest_point()
    fieldMesh_x,fieldMesh_y = np.meshgrid(lons,lats)
    landMask = np.loadtxt('datafile_landMask_677x_253y')
    coastMask = np.loadtxt('datafile_coastMask_677x_253y')
    
    landBorder = landMask.copy()
    val_add = 2
    for i1 in range(6):
        landBorder = getLandBorder(landBorder,lons,lats,val_add)
        val_add += 1
    
    def closest_index(lat,lon,mask_test):
        distMat = 1e5 * np.ones(fieldMesh_x.shape)
        
        test_indices = np.where(mask_test == 1)
        
        distMat_lon = (lon - fieldMesh_x[test_indices[0],test_indices[1]])*1.11e2*0.788 #find distances coastal element w.r.t. ocean cells. 0.788 comes from lat=38deg (medit.)
        distMat_lat = (lat - fieldMesh_y[test_indices[0],test_indices[1]])*1.11e2
    
        distMat[test_indices[0],test_indices[1]] = np.sqrt(np.power(distMat_lon, 2) + np.power(distMat_lat, 2))    
        
        return np.where(distMat == distMat.min())[0][0],np.where(distMat == distMat.min())[1][0]
    
    ### interpolate beaching to closest coastal cell
    hist_beaching_coast = np.zeros(fieldMesh_x.shape)
    for i1 in range(len(lats)):
        for i2 in range(len(lons)):
            
            if hist_beaching_sum[i1,i2] > 0:
                
                i_lat,i_lon = closest_index(lats[i1],lons[i2],coastMask)
                hist_beaching_coast[i_lat,i_lon] += hist_beaching_sum[i1,i2]
    
    
    hist_beaching_coast /= grid_area
    
    #%%
    ### go through the landborder defined above with increased width
    hist_beaching_extended = np.zeros(fieldMesh_x.shape)
    indices_border = np.where(landBorder > 1)            
    for i1 in range(len(indices_border[0])):
        lon_ = lons[indices_border[1][i1]]
        lat_ = lats[indices_border[0][i1]]
        i_lat,i_lon = closest_index(lat_,lon_,coastMask)
        
        hist_beaching_extended[indices_border[0][i1],indices_border[1][i1]] += hist_beaching_coast[i_lat,i_lon]
    
    #%%    
    hist_sinking_area = hist_sinking_sum / grid_area
    hist_flux_beaching = hist_beaching_extended/len(date_arr)
    hist_flux_sinking = hist_sinking_area/len(date_arr)
    
    cmap1 = cmap=plt.cm.viridis_r
    
    cmap2 = cmap=cmocean.cm.turbid
    
    fig = plt.figure(figsize=(7,5),dpi=120)     
    ax = plt.axes(projection=ccrs.PlateCarree())
    sea_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                            edgecolor='grey',
                                            facecolor='silver',zorder=0)
    ax.add_feature(sea_50m)


    log_min_2 = -5
    log_max_2 = 0.8*np.log10(hist_flux_sinking).max()
    cplot2 = ax.contourf(fieldMesh_x,fieldMesh_y,hist_flux_sinking,levels=np.logspace(-5,-1,9),cmap=cmap2,norm=colors.LogNorm(),transform=ccrs.PlateCarree(),extend='both')

    log_min_1 = -3
    log_max_1 = np.log10(hist_flux_beaching).max()
    cplot = ax.contourf(fieldMesh_x,fieldMesh_y,hist_flux_beaching,levels=np.logspace(-3,1,9),cmap=cmap1,norm=colors.LogNorm(),transform=ccrs.PlateCarree(),extend='both')
        
    cax1 = fig.add_axes([0.15, 0.25, 0.7, 0.02])
    cax2 = fig.add_axes([0.15, 0.13, 0.7, 0.02])
    
    cbar1 = plt.colorbar(cplot,cax=cax1,orientation='horizontal',ticks=[10**i for i in np.linspace(log_min_1,int(log_max_1),int(log_max_1)-log_min_1+1)])
    cbar2 = plt.colorbar(cplot2,cax=cax2,orientation='horizontal',ticks=[10**i for i in np.linspace(log_min_2,int(log_max_2),int(log_max_2)-log_min_2+1)])

    cbar1.ax.tick_params(labelsize=7)
    cbar1.set_label('beaching, 2006-2016 [kg/km$^2$/day]',fontsize=7)
    cbar2.ax.tick_params(labelsize=7)
    cbar2.set_label('sinking, 2006-2016 [kg/km$^2$/day]',fontsize=7)
    
    
               
            