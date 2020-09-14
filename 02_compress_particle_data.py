#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 16:23:11 2018
Parcels gives a set of netcdf files containing all particle locations at a given dt.
Not all data is necessary: not each day contains a measurement. This script combines all necessary particle data for a given measurement: 
it combines particles from different sources (population,rivers, fisheries) and stores all relevant information in .csv files

1)Read plastic measurements from the plastic_data excel sheet
2)Read particle files as created by parcels
3)Look up measurement times in particle files
4)Files are written for each measurement for further analysis/optimization

@author: kaandorp
"""

from datetime import timedelta
from datetime import datetime
import numpy as np
import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
import os
import xarray as xr
import time


def calculateWindSpeed(t_curr,lon,lat):
    """
    Accurate to 3hrs    
    """
    yr = '%i' % t_curr.year
    mon = '%2.2i' % t_curr.month 

    file_str = '*%s%s_wnd.nc' % (yr,mon)
    filename = glob.glob(os.path.join('/Users/kaandorp/Data/Wind/WW3',yr,file_str))[0]
    
    data_w = Dataset(filename)
    
    data_time = data_w['time'][:]
    data_date = np.array([datetime(1990,1,1) + timedelta(days = dt) for i,dt in enumerate(data_time.data[:])])
    
    index_closest = np.argmin(np.abs(t_curr - data_date))
 
    X,Y = np.meshgrid(data_w['longitude'][:].data,data_w['latitude'][:].data)
    distMat = (X-lon)**2 + (Y-lat)**2
    distMat = np.ma.array(distMat,mask=data_w['uwnd'][0,:,:].mask)
    
    i_min_lon = np.where(distMat == distMat.min())[1][0]
    i_min_lat = np.where(distMat == distMat.min())[0][0]
    
    w_closest = np.sqrt((data_w['uwnd'][index_closest,i_min_lat,i_min_lon])**2 + (data_w['vwnd'][index_closest,i_min_lat,i_min_lon])**2)
    
    if np.isnan(w_closest):
        raise RuntimeError('!!!')
    
    return w_closest




plt.close('all')



#%% Init

units_iter = ['g/km2','#/km2']

Kval = 10

dataDir = 'noDelete/K%i' %Kval #cartesius

projectFolder = '34_CSVFiles_noDelete_popRivFis_K%i/' %Kval #cartesius


use_external_SSD = True
save_particle_files = True #files which contain all particle information for a given day
save_measurement_files = False #files which contain measurement data for a given day

use_pop = True #the different sources
use_rivers = True
use_fisheries = True

n_folders_per_type = 6 #amount of folders used per type, can be set lower for less data and more speed



if os.environ['USER'] == 'kaandorp': # desktop
    if use_external_SSD:
        homeDir = '/Volumes/externe_SSD/kaandorp'
    else:
        homeDir = os.environ['HOME']
    
    outDir = os.path.join(homeDir, projectFolder)
    dataFile = os.path.join(homeDir, 'Documents', 'plastic_data_v4.xlsx')
    
    if use_pop:
        particleFolders_pop = os.path.join(homeDir , dataDir, 'popOnly/')
        folders_pop = np.sort(glob.glob(os.path.join(particleFolders_pop,'*')))[0:n_folders_per_type]
    if use_rivers:
        particleFolders_rivers = os.path.join(homeDir , dataDir, 'riversOnly/')
        folders_rivers = np.sort(glob.glob(os.path.join(particleFolders_rivers,'*')))[0:n_folders_per_type]
    if use_fisheries:
        particleFolders_fisheries = os.path.join(homeDir , dataDir, 'fisheriesOnly/')
        folders_fisheries = np.sort(glob.glob(os.path.join(particleFolders_fisheries,'*')))[0:n_folders_per_type]

elif os.environ['USER'] == 'mikaelk': # cartesius
    homeDir = '/scratch-shared/mikaelk'
  
    outDir = os.path.join(homeDir, projectFolder)
    dataFile = os.path.join(homeDir, 'Documents', 'plastic_data_v4.xlsx')
    
    if use_pop:
        particleFolders_pop = os.path.join(homeDir , dataDir, 'popOnly/')
        folders_pop = np.sort(glob.glob(os.path.join(particleFolders_pop,'*')))[0:n_folders_per_type]
    if use_rivers:
        particleFolders_rivers = os.path.join(homeDir , dataDir, 'riversOnly/')
        folders_rivers = np.sort(glob.glob(os.path.join(particleFolders_rivers,'*')))[0:n_folders_per_type] 
    if use_fisheries:
        particleFolders_fisheries = os.path.join(homeDir , dataDir, 'fisheriesOnly/')
        folders_fisheries = np.sort(glob.glob(os.path.join(particleFolders_fisheries,'*')))    [0:n_folders_per_type]
    
else:
    raise RuntimeError('unknown user environment')
    


if os.path.exists(outDir):
    print('Writing files to %s\n' % outDir)
else:
    os.makedirs(outDir)
    print('WARNING: Creating folder %s\n' % outDir)




#------------- file containing measurements -------------------
data = pd.read_excel(dataFile,
sheet_name=0,
header=0,
index_col=False,
keep_default_na=True,
dtype={'depth':str,'method':str}
)

#-------------Cleaning and selecting measurement data-----------------------

#removing uncertainties
data['size_class'] = data['size_class'].str.replace('?','')
data['Kukulka_corr'] = data['Kukulka_corr'].str.replace('?','')


data_d = data['depth'].str.replace(',','.')
is_unknown_d = (data_d == '-')
is_estim_d = data_d.str.contains('\?',na=False)
is_known_d = ~is_unknown_d & ~is_estim_d
data_d =data_d.str.replace('?','')

data.loc[:,'depth'] = pd.to_numeric(data_d,errors='coerce')


is_manta = data['method'].str.contains('manta|neuston')
trawl_sizes = data['method'].str.extract('.*\((.*)\).*')

trawl_height    = data['method'].str.split('(').str[1].str.split(')').str[0].str.split('x').str[0].astype(np.float64)
trawl_width     = data['method'].str.split('(').str[1].str.split(')').str[0].str.split('x').str[1].astype(np.float64)


is_unknown_h = np.isnan(trawl_height)

bool_i = is_manta & is_unknown_h & is_unknown_d
data.loc[bool_i,'depth'] = 0.1
bool_i = is_manta & is_unknown_d & (trawl_height <= 0.2)
data.loc[bool_i,'depth'] = 0.1
bool_i = is_manta & is_unknown_d & (trawl_height > 0.2)
data.loc[bool_i,'depth'] = 0.15


for selected_unit in units_iter:
    
    if selected_unit == 'g/km2':
        string_meastype = 'weight'
    elif selected_unit == '#/km2':
        string_meastype = 'num'
    else:
        raise RuntimeError('unknown unit type')
    

    data_selectedUnit = data[( (data['size_class'] == '<5mm')  | (data['size_class'] == 'all') ) & (~np.isnan(data[selected_unit])) ]
    
    
    #%% Main program
    
    t_start = time.time()
    
    data_analysis_check = {}
    data_analysis_check['KDE'] = np.array([])
    
    data_analysis={}
    data_analysis['KDE'] = np.array([])
    data_analysis['measured'] = np.array([])
    data_analysis['studyID'] = np.array([])
    data_analysis['lon'] = np.array([])
    data_analysis['lat'] = np.array([])
    
    list_concat = []
    #TODO: check whether all folders have all time steps, and what amount of files
    if use_pop:
        folder_0_pop = folders_pop[0]
        fileList_p0_pop = list(np.sort(glob.glob(os.path.join(folder_0_pop ,'K%3.1f_data*' % Kval))))
        list_concat.append(folders_pop)
    if use_rivers:
        folder_0_rivers = folders_rivers[0] 
        fileList_p0_rivers = list(np.sort(glob.glob(os.path.join(folder_0_rivers ,'K%3.1f_data*' % Kval))))
        list_concat.append(folders_rivers)
    if use_fisheries:
        folder_0_fisheries = folders_fisheries[0]
        fileList_p0_fisheries = list(np.sort(glob.glob(os.path.join(folder_0_fisheries ,'K%3.1f_data*' % Kval))))
        list_concat.append(folders_fisheries)
    
    n_particleFiles = len(fileList_p0_pop)
    
    
    files_initStates = [glob.glob(os.path.join(folder,'particlesRelease*')) for folder in np.concatenate(list_concat)]
    
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
    init_lat = initStates.iloc[:,1]
    init_lon = initStates.iloc[:,2]
    
    
    fileCounter = 0
    
    # verify correct implementation for looking up particle origin
    lat_origin = []
    lon_origin = []
    lat_ptcl = []
    lon_ptcl = []
    
    
    
    for i1 in range(100,n_particleFiles): #mass: start at 138    #: start at 115
        list_stack = []
        list_concat = []
        
        if use_pop:
            #--------------------- open particles due to population density-----------------------
            fileName_p = os.path.basename(fileList_p0_pop[i1])
            pFiles = [os.path.join(folders_pop[i2],fileName_p) for i2 in range(len(folders_pop))]
            
            # look if files exist (due to missing files from parcels memory allocation error)
            fileExist = [os.path.isfile(pFiles[i2]) for i2 in range(len(pFiles))]
            pFiles = list(np.array(pFiles)[fileExist])
            
            data_particles_pop = xr.open_mfdataset(pFiles,concat_dim='traj',decode_times=False)
        
            traj = data_particles_pop['trajectory']
            locs_zero = np.array([0])
            locs_zero = np.append(locs_zero,np.where((traj[2:,0] > traj[1:-1,0]) & (traj[0:-2,0] > traj[1:-1,0]))[0] + 1)
         
            index_start = 0
            index_end = len(folders_pop)
            p_indices_p = data_particles_pop.trajectory[:,0].compute()
            for i2 in range(len(pFiles)):
                if i2 < len(pFiles)-1:
                    p_indices_p[locs_zero[i2]:locs_zero[i2+1]] = p_indices_p[locs_zero[i2]:locs_zero[i2+1]] + index_plus[index_start:index_end][fileExist][i2] #*index_plus #i2*index_plus
                else:
                    p_indices_p[locs_zero[i2]:] = p_indices_p[locs_zero[i2]:] + index_plus[index_start:index_end][fileExist][i2] #np.arange(len(fileExist))[fileExist][i2]*index_plus
            p_indices_p = p_indices_p.values
            list_stack.append(p_indices_p)
            list_concat.append(data_particles_pop)
        
        if use_rivers:
            #-------------------- open particles due to rivers -----------------------------------
            fileName_p = os.path.basename(fileList_p0_rivers[i1])
            pFiles = [os.path.join(folders_rivers[i2],fileName_p) for i2 in range(len(folders_rivers))]
            
            # look if files exist (due to missing files from parcels memory allocation error)
            fileExist = [os.path.isfile(pFiles[i2]) for i2 in range(len(pFiles))]
            pFiles = list(np.array(pFiles)[fileExist])
            
            data_particles_rivers = xr.open_mfdataset(pFiles,concat_dim='traj',decode_times=False)
        
            traj = data_particles_rivers['trajectory']
            locs_zero = np.array([0])
            locs_zero = np.append(locs_zero,np.where((traj[2:,0] > traj[1:-1,0]) & (traj[0:-2,0] > traj[1:-1,0]))[0] + 1)
         
            index_start = len(folders_pop)
            index_end = len(folders_pop)+len(folders_rivers)
            p_indices_r = data_particles_rivers.trajectory[:,0].compute()
            for i2 in range(len(pFiles)):
                if i2 < len(pFiles)-1:
                    p_indices_r[locs_zero[i2]:locs_zero[i2+1]] = p_indices_r[locs_zero[i2]:locs_zero[i2+1]] + index_plus[index_start:index_end][fileExist][i2] #*index_plus #i2*index_plus
                else:
                    p_indices_r[locs_zero[i2]:] = p_indices_r[locs_zero[i2]:] + index_plus[index_start:index_end][fileExist][i2] #np.arange(len(fileExist))[fileExist][i2]*index_plus
            p_indices_r = p_indices_r.values
            list_stack.append(p_indices_r)
            list_concat.append(data_particles_rivers)

        if use_fisheries:
            #----------------------open particles due to fisheries -------------------------------
            fileName_p = os.path.basename(fileList_p0_fisheries[i1])
            pFiles = [os.path.join(folders_fisheries[i2],fileName_p) for i2 in range(len(folders_fisheries))]
            
            # look if files exist (due to missing files from parcels memory allocation error)
            fileExist = [os.path.isfile(pFiles[i2]) for i2 in range(len(pFiles))]
            pFiles = list(np.array(pFiles)[fileExist])
            
            data_particles_fisheries = xr.open_mfdataset(pFiles,concat_dim='traj',decode_times=False)
            
            traj = data_particles_fisheries['trajectory']
            locs_zero = np.array([0])
            locs_zero = np.append(locs_zero,np.where((traj[2:,0] > traj[1:-1,0]) & (traj[0:-2,0] > traj[1:-1,0]))[0] + 1)
         
            index_start = len(folders_pop)+len(folders_rivers)
            index_end = -1
            p_indices_f = data_particles_fisheries.trajectory[:,0].compute()
            for i2 in range(len(pFiles)):
                if i2 < len(pFiles)-1:
                    p_indices_f[locs_zero[i2]:locs_zero[i2+1]] = p_indices_f[locs_zero[i2]:locs_zero[i2+1]] + index_plus[index_start:index_end][fileExist][i2] #*index_plus #i2*index_plus
                else:
                    p_indices_f[locs_zero[i2]:] = p_indices_f[locs_zero[i2]:] + index_plus[index_start:index_end][fileExist][i2] #np.arange(len(fileExist))[fileExist][i2]*index_plus
            p_indices_f = p_indices_f.values
            list_stack.append(p_indices_f)
            list_concat.append(data_particles_fisheries)
        
        p_indices = np.hstack(list_stack)
    
    
    
        #--------------------- combine data --------------------------------------------------
        data_particles = xr.concat(list_concat,dim='traj')
        
    
        if data_particles.time.units == 'seconds since 2001-01-01T00:00:00.000000000':
            #compute which particle contains all times, use this one to create the time array
            index_valid = np.where(~np.isnan(data_particles.time.data[:,data_particles.time.data.shape[1]-1]))[0][-1]
            t_range = np.sort(np.array([datetime(2001,1,1) + timedelta(seconds = int(data_particles.time.data[index_valid,i1].compute())) for i1 in range(len(data_particles.time.data[index_valid,:])) ]))
    #        t_range = pd.to_datetime(np.sort(data_particles.time.data[0,:]))
            if (t_range[-1] - t_range[-2]).days != 1:
                print(t_range)
        
        else:
           raise RuntimeError('check time units')
           
           
        # check which measurements are relevant
        selectedMeasurements = (data_selectedUnit['t'] > t_range[0]) & (data_selectedUnit['t'] <= t_range[-1])
        
        if selectedMeasurements.any():
            print('File %s contains useful information' % fileList_p0_pop[i1])
            
            data_measurements_inTimeRange = data_selectedUnit[selectedMeasurements]
            # Sort data_measurements_inTimeRange by time
            data_measurements_inTimeRange = data_measurements_inTimeRange.sort_values(by=['t'])
            

            for index_meas in data_measurements_inTimeRange.index:    
                
                t_curr = data_measurements_inTimeRange['t'][index_meas]
    
                # select closest time (or two closest times)
                t_delta = t_curr - t_range
                index_closest = np.argmin(np.abs(t_delta))
                t_closest = t_range[index_closest]
                

                a = np.where(np.isnan(data_particles.lon.data[:,index_closest]))
                b = np.where(np.isnan(data_particles.lat.data[:,index_closest]))
                c = np.where(np.isnan(data_particles.time.data[:,index_closest]))
                indices_delete = np.unique(np.concatenate((a,b,c),0))
                
                particles_valid={}
                particles_valid['ID'] = np.delete(p_indices,indices_delete)
                particles_valid['lon'] = np.delete(data_particles.lon.data[:,index_closest],indices_delete)
                particles_valid['lat'] = np.delete(data_particles.lat.data[:,index_closest],indices_delete)
                
                particles_valid['coastalAge_1cell'] = np.delete(data_particles.coastalAge_1cell.data[:,index_closest],indices_delete)
                particles_valid['coastalAge_2cell'] = np.delete(data_particles.coastalAge_2cell.data[:,index_closest],indices_delete)
                particles_valid['age'] = np.delete(data_particles.age.data[:,index_closest],indices_delete)
                
                particles_valid['dot_curr'] = np.delete(data_particles.dot_curr.data[:,index_closest],indices_delete)
                particles_valid['dot_wind'] = np.delete(data_particles.dot_wind.data[:,index_closest],indices_delete)
                particles_valid['dot_wind_Ek'] = np.delete(data_particles.dot_wind_Ek.data[:,index_closest],indices_delete)
                particles_valid['dot_stokes'] = np.delete(data_particles.dot_stokes.data[:,index_closest],indices_delete)
    
                particles_valid['mag_curr'] = np.delete(data_particles.mag_curr.data[:,index_closest],indices_delete)
                particles_valid['mag_stokes'] = np.delete(data_particles.mag_stokes.data[:,index_closest],indices_delete)
                particles_valid['mag_wind'] = np.delete(data_particles.mag_wind.data[:,index_closest],indices_delete)
    

                if save_measurement_files:
                    U_10 = calculateWindSpeed(t_curr,data_measurements_inTimeRange['lon'][index_meas],data_measurements_inTimeRange['lat'][index_meas])
                    
                    U_10_bft = (U_10/0.837)**(2/3)
                    wave_frombft = (0.3/9.81)*U_10_bft**2 # relation from Rossby and Montgomery (1935)
                    
                    wb = 5.3e-3
                    U10_to_fric = 1.2e-3
                    k = 0.4
                    
                    def Kukulka_corr(U10,VHM0,depth):
                        u_fric = U10_to_fric*U10
                        
                        corr = u_fric > 0.006
                        
                        if corr:
                            return_val = 1/(1 - np.exp(-depth * wb * (1/(1.5*u_fric*k*VHM0)) )) 
                        else:
                            return_val = 1.0
                        
                        return return_val
        
                    if data_measurements_inTimeRange['Kukulka_corr'][index_meas] == 'n':
                        Kukulka_factor = Kukulka_corr(U_10,wave_frombft,data_measurements_inTimeRange['depth'][index_meas])
                    else:
                        Kukulka_factor = 1.0
                    if np.isnan(Kukulka_factor):
                        print('nan')
                        print(i1)
                        print(index_meas)
                    
                    #Size range correction factor ('<5mm' vs 'all' categories for g/km2)
                    sizeRange_correction = 1.
                    if data_measurements_inTimeRange['size_class'][index_meas] == '<5mm':
                        
                        if selected_unit == 'g/km2':        
                            
                            print('Weight correction added, %s' % data_measurements_inTimeRange['key'][index_meas])
                            sizeRange_correction = 15.46
                        if selected_unit == '#/km2':
                            sizeRange_correction = 1.14
        
        
                    dataArray_measurement = np.array([data_measurements_inTimeRange['study_ID'][index_meas],data_measurements_inTimeRange['lon'][index_meas],data_measurements_inTimeRange['lat'][index_meas],
                                                      (t_curr-datetime(2001,1,1)).total_seconds(),data_measurements_inTimeRange[selected_unit][index_meas],Kukulka_factor,sizeRange_correction])
                    dataArray_measurement_header = ','.join(['study_ID','lon','lat','dt_sec(2001-1-1)','measured density','Kukulka factor','size range factor'])
                    
                    c_meas = 0
                    new_filename_required = True
                    while new_filename_required:
                        measurement_filename = os.path.join(outDir,'%i%2.2i%2.2i_measurement_%s_%3.3i.csv' % (t_closest.year,t_closest.month,t_closest.day,string_meastype,c_meas) )
                        if os.path.exists(measurement_filename):
                            c_meas += 1
                        else:
                            np.savetxt(measurement_filename, dataArray_measurement[None],delimiter=',',newline='\n',header=dataArray_measurement_header)
                            new_filename_required = False
                        
                if save_particle_files:       
                    dataArray_particles = np.array([particles_valid[key] for key in particles_valid.keys()])
                    dataArray_particles_header = ','.join(particles_valid.keys())
                    
                    particles_filename = os.path.join(outDir,'%i%2.2i%2.2i_particles.csv' % (t_closest.year,t_closest.month,t_closest.day) )
                    if os.path.exists(particles_filename):
                        pass
                    else:
                        np.savetxt(particles_filename,dataArray_particles.T,delimiter=',',newline='\n',header=dataArray_particles_header)
    
                fileCounter += 1
    
    
    t_total = time.time() - t_start
    print(t_total)
    
                      