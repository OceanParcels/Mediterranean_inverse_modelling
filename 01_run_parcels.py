#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 07 16:19:35 2018
Script to run parcels. in the command line, use for example 'python 01_run_parcels.py -k 10 -f output_folder -p 0.6 -w 0.2 -n 5'
@author: kaandorp
"""


import sys, getopt
import shapefile
import csv
from parcels import FieldSet, VectorField, ParticleSet, JITParticle, ErrorCode, Field, Variable
from datetime import timedelta
from datetime import datetime
import calendar
import numpy as np
import glob
from netCDF4 import Dataset
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from parcels import rng as random
import math
import pandas as pd
import cartopy.crs as ccrs
import os
import matplotlib
from parcels.tools.converters import Geographic, GeographicPolar 
import time

plt.switch_backend('agg')

def getArguments(argv):
    """
    Function to parse command line arguments
    """
    
    K = 10.0
    fname = '01_defaultFolder'
    fraction_pop = 0.6
    fraction_fisheries = 0.2
    particles_per_day = 5

    try:
        opts, args = getopt.getopt(argv,"hk:f:p:w:n:")

    except getopt.GetoptError:
        print('main.py -k <viscosity value>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Usage:\n main.py -k <viscosity value> -f <outputfolder name> -p <waste due to pop. density> ')
        elif opt == '-k':
            K = float(arg)
        elif opt == '-f':
            fname = str(arg)
        elif opt == '-p':
            fraction_pop = float(arg)            
        elif opt == '-w':
            fraction_fisheries = float(arg)
        elif opt == '-n':
            particles_per_day = int(arg)
            
    return K,fname,fraction_pop,fraction_fisheries,particles_per_day


def getLandMask(filename, fieldname):
    """
    Get mask of the land using a netcdf file.
    filename: name of netcdf file
    fieldname: field from which to obtain land mask
    """
    cdffile = Dataset(filename, 'r')
    f = cdffile.variables[fieldname][:]
    
    L= f[0,0,:,:].mask    
    L = L*1
    return L


def getCoastMask(landMask,lon,lat):
    """
    Function to obtain a mask of the coast, uses the landmask and searches for 
    boundaries with the sea (horiz./vert. adjacent cells only)
    TODO: check for cyclic boundaries
    TODO: check diagonal cells as well?
    """
    
    n_lat = landMask.shape[0]
    n_lon = landMask.shape[1]
    
    coastMask = np.zeros([n_lat,n_lon])

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
                    if landMask[i1-1,i2] == 0:
                        coastMask[i1-1,i2] = 1
                if check_bot:
                    if landMask[i1+1,i2] == 0:
                        coastMask[i1+1,i2] = 1
                if check_left:
                    if landMask[i1,i2-1] == 0:
                        coastMask[i1,i2-1] = 1
                if check_right:
                    if landMask[i1,i2+1] == 0:
                        coastMask[i1,i2+1] = 1
            
    return coastMask


def getLandBorder(landMask,lon,lat):
    """
    Function to obtain a mask of the land which borders ocrean, uses the landmask and searches for 
    boundaries with the sea (horiz./vert. adjacent cells only)
    TODO: check for cyclic boundaries
    TODO: check diagonal cells as well?
    """
    
    n_lat = landMask.shape[0]
    n_lon = landMask.shape[1]
    
    borderMask = np.zeros([n_lat,n_lon])

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
                    if landMask[i1-1,i2] == 0:
                        borderMask[i1,i2] = 1
                if check_bot:
                    if landMask[i1+1,i2] == 0:
                        borderMask[i1,i2] = 1
                if check_left:
                    if landMask[i1,i2-1] == 0:
                        borderMask[i1,i2] = 1
                if check_right:
                    if landMask[i1,i2+1] == 0:
                        borderMask[i1,i2] = 1
            
    return borderMask


def getParticleCoordinates(pset,time_sec):
    lonArr = np.array([pset[i1].lon for i1 in range(len(pset))])
    latArr = np.array([pset[i1].lat for i1 in range(len(pset))])
    timeArr = np.array([pset[i1].time for i1 in range(len(pset))])
    
    lonArr = lonArr[timeArr <= time_sec]
    latArr = latArr[timeArr <= time_sec]
    return np.array([lonArr,latArr])
    

def plotParticles(pset,fieldMesh_x,fieldMesh_y,plottingDomain,landMask,time_sec,figsize=(30,12),res_x=200,res_y=150):
    

    plt.figure(figsize=figsize)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_extent((plottingDomain[0],plottingDomain[1],plottingDomain[2],plottingDomain[3]), ccrs.PlateCarree())
    partCoordinates = getParticleCoordinates(pset,time_sec)
    ax.plot(partCoordinates[0,:],partCoordinates[1,:],'ko',markersize=5,transform=ccrs.PlateCarree())
    
#    res_x = 200
#    res_y = 150
    
    #plot first set of the fieldlist
    if len(pset.fieldset.UV.U.data) > 1:
        U = pset.fieldset.UV.U.data[0]
        V = pset.fieldset.UV.V.data[0]
    else:
        U = pset.fieldset.UV.U.data
        V = pset.fieldset.UV.V.data
        
    plotGrid_x, plotGrid_y = np.meshgrid(np.linspace(fieldMesh_x.min(),fieldMesh_x.max(),res_x),np.linspace(fieldMesh_y.min(),fieldMesh_y.max(),res_y))
    U_plot = griddata((fieldMesh_x.ravel(),fieldMesh_y.ravel()),U[:,:].ravel(),(plotGrid_x,plotGrid_y))
    V_plot = griddata((fieldMesh_x.ravel(),fieldMesh_y.ravel()),V[:,:].ravel(),(plotGrid_x,plotGrid_y))
    landMask_plot = griddata((fieldMesh_x.ravel(),fieldMesh_y.ravel()), landMask.ravel(),(plotGrid_x,plotGrid_y))
    speed = np.sqrt(U_plot**2+V_plot**2)
    
    U_masked = np.ma.array(U_plot,mask=landMask_plot.astype(bool))
    V_masked = np.ma.array(V_plot,mask=landMask_plot.astype(bool))
    #ax.quiver(plotGrid_x, plotGrid_y,U_masked,V_masked,speed,cmap=plt.cm.gist_ncar,scale=10,transform=ccrs.PlateCarree())
    ax.quiver(plotGrid_x, plotGrid_y,U_masked,V_masked,speed,cmap=plt.cm.gist_ncar,scale=20,width=0.001,transform=ccrs.PlateCarree())


def initializePopulationMatrices_50km(pop_dataFile,coastMask,fieldMesh_x,fieldMesh_y,MPWFile='./datafile_MPW_gridded'):
    """
    Initialize the coastal population density matrices. This is one matrix for every specified year. 
    MPW density in a 50km range is considered as done in e.g. Lebreton 2018
    These matrices can then later be linearly interpolated for densities at a given time
    """
    mpw_gridded = np.loadtxt(MPWFile)

    populationCoastMatrix = np.zeros([5,fieldMesh_x.shape[0],fieldMesh_x.shape[1]])
    pop_datafile            = os.path.join(homeDir , 'Data/populationDensity/gpw_v4_e_atotpopbt_dens_2pt5_min.nc')
    
    dataPop = Dataset(pop_datafile)
    ids_long = [4175,5208]
    ids_lat = [935,1320]
    
    for i1 in range(5):
        meshPop_x, meshPop_y = np.meshgrid(dataPop.variables['longitude'][ids_long[0]:ids_long[1]],dataPop.variables['latitude'][ids_lat[0]:ids_lat[1]])
        densityPop = dataPop.variables['Population Density, v4.10 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'][i1,ids_lat[0]:ids_lat[1],ids_long[0]:ids_long[1]]
        
    
        densityPop_i = griddata((meshPop_x.ravel(),meshPop_y.ravel()),densityPop.ravel(),(fieldMesh_x,fieldMesh_y))
        densityPop_i[densityPop_i < 0] = 0
        densityPop_i[np.isnan(densityPop_i)] = 0 
    
        
        densityPop_i = densityPop_i * mpw_gridded
        
        indices_coast = np.where(coastMask == 1)
        ids_lon_ = indices_coast[1]
        ids_lat_ = indices_coast[0]
        
        for i2 in range(len(indices_coast[0])):
            lon_ = lons[ids_lon_[i2]]
            lat_ = lats[ids_lat_[i2]]      
            
            id_lon_min = max((ids_lon_[i2]-8),0)
            id_lat_min = max((ids_lat_[i2]-8),0)
            
            id_lon_max = min((ids_lon_[i2]+9),677)
            id_lat_max = min((ids_lat_[i2]+9),253)
            
            mesh_x_centered = fieldMesh_x[id_lat_min:id_lat_max,id_lon_min:id_lon_max]
            mesh_y_centered = fieldMesh_y[id_lat_min:id_lat_max,id_lon_min:id_lon_max]
            densityPop_centered = densityPop_i[id_lat_min:id_lat_max,id_lon_min:id_lon_max]
            
            
            distance_mat = np.sqrt( ((lon_ - mesh_x_centered)*np.cos(lat_*np.pi/180)*1.11e2)**2 + ((lat_ - mesh_y_centered)*1.11e2)**2)
            id_dist_50 = (distance_mat < 50)
            
            sum_50km = densityPop_centered[id_dist_50].sum()
            
            populationCoastMatrix[i1,ids_lat_[i2],ids_lon_[i2]] += sum_50km
    
    return populationCoastMatrix
        
    
def initializeRiverMatrices(riverShapeFile, pollutionFile, coastMask, fieldMesh_x, fieldMesh_y, selection_estim='mid', plot_coastalinput = False):
    """
    Initialize the river release matrices. This is one matrix for every specifed month.
    Based on data from Lebreton, where catchment areas/runoff were estimated together with
    MPW from Jambeck
    selection_estim: choices are 'mid','low','high', which correspond to the 
    lower-upper and middle of the confidence bounds of MPW emitted by rivers from
    Lebreton
    """

    print('Initializing river input matrices...')

    #import shapefile
    sf = shapefile.Reader(riverShapeFile)
    
    #extract files within mediterranean
    plottingDomain = [-6,37,30,46]
    
    rivers = {}
    rivers['longitude'] = np.array([])
    rivers['latitude'] = np.array([])
    rivers['ID'] = np.array([],dtype=int)
    rivers['dataArray'] = np.array([])
    
    for i1 in range(len(sf.shapes())):
        long = sf.shape(i1).points[0][0]
        lat = sf.shape(i1).points[0][1]
        
        if plottingDomain[0] < long <plottingDomain[1] and plottingDomain[2] < lat < plottingDomain[3]:
            rivers['longitude'] = np.append(rivers['longitude'],long)
            rivers['latitude'] = np.append(rivers['latitude'],lat)
            rivers['ID'] = np.append(rivers['ID'],i1)
            
            
    with open(pollutionFile, 'r',encoding='ascii') as csvfile:
        filereader = csv.reader(csvfile, delimiter=';')
        i1 = 0
        for row in filereader:
            
            if i1 == 0:
                riverHeaders = row
            
            if i1 > 0:
                
            
                data_ID = i1-1
                
                if i1 == 1:
                    dataArray = [float(row[i2].replace(',','.')) for i2 in range(len(row))]
                    rivers['dataArray'] = dataArray
                else:
                    if data_ID in rivers['ID']:
                        dataArray = [float(row[i2].replace(',','.')) for i2 in range(len(row))]
                        rivers['dataArray'] = np.vstack([rivers['dataArray'],dataArray])
            i1 += 1
        
    # check which columns contain data for the specified confidence level from selection_estim parameter    
    columnNumbers = []
#    i1 = 0
    for idx, columnName in enumerate(riverHeaders):
        if (selection_estim + '_') in columnName:
            columnNumbers.append(idx)
#        i1 += 1
    assert(len(columnNumbers) == 12), "there should be 12 entries corresponding to waste emitted per month"
    
    coastIndices = np.where(coastMask == 1)
    assert(np.shape(coastIndices)[0] == 2), "coastMask.data should be 2 by something"
    
    # array containing indices of rivers not belonging to mediterranean, which are to be deleted
    deleteEntries = np.array([],dtype=int)
    
    # matrix corresponding to fieldmesh, with per coastal cell the amount of river pollution
    riverInputMatrix = np.zeros([12,fieldMesh_x.shape[0],fieldMesh_x.shape[1]])
    
    # for every river
    for i1 in range(len(rivers['longitude'])):
        
        lon_river = rivers['longitude'][i1]
        lat_river = rivers['latitude'][i1]
        
        dist = 1e10
        # check which point is closest
        for i2 in range(np.shape(coastIndices)[1]):
            lon_coast = lons[coastIndices[1][i2]]
            lat_coast = lats[coastIndices[0][i2]]
        
            lat_dist = (lat_river - lat_coast) * 1.11e2
            lon_dist = (lon_river - lon_coast) * 1.11e2 * np.cos(lat_river * np.pi / 180)
            dist_tmp = np.sqrt(np.power(lon_dist, 2) + np.power(lat_dist, 2))
            
            # save closest distance
            if dist_tmp < dist:
                dist = dist_tmp
                lat_ID = coastIndices[0][i2]
                lon_ID = coastIndices[1][i2]
            
        # if distance to closest point > threshold (3*approx cell length), delete entry
        if dist > 3*0.125*1.11e2:
            deleteEntries = np.append(deleteEntries,i1)
        # else: get pollution river, and add to releasematrix
        else:
            # add plastic input as obtained from the dataset
            for idx, val in enumerate(columnNumbers):
                riverInputMatrix[idx,lat_ID,lon_ID] += rivers['dataArray'][i1,val]
            
    
    # rivers ending in mediterranean
    rivers_medit = {}
    rivers_medit['longitude'] = np.delete(rivers['longitude'],deleteEntries)
    rivers_medit['latitude'] = np.delete(rivers['latitude'],deleteEntries)
    rivers_medit['ID'] = np.delete(rivers['ID'],deleteEntries)
    rivers_medit['dataArray'] = np.delete(rivers['dataArray'],deleteEntries,axis=0)
            
    if plot_coastalinput:
        times = ['jan','feb','mar','apr','may','june','july','aug','sep','oct','nov','dec']
        for i1 in range(riverInputMatrix.shape[0]):
            
            minval = 1e-15
            maxval = riverInputMatrix.max()
            riverInputMatrix_plt = np.copy(riverInputMatrix)
            riverInputMatrix_plt[riverInputMatrix_plt == 0.0] = 1e-15
            figsize=(30,12)
            plt.figure(figsize=figsize)
            ax = plt.axes(projection=ccrs.PlateCarree())
            ax.coastlines(resolution='10m')
            cmap = plt.cm.inferno
            handle = ax.contourf(fieldMesh_x,fieldMesh_y,riverInputMatrix_plt[i1,:,:],vmin=minval,vmax=maxval,norm=matplotlib.colors.LogNorm(),cmap=cmap,transform=ccrs.PlateCarree())
            plt.colorbar(handle)
            plt.title('River inputs, %s' % times[i1])
    
    return riverInputMatrix

def populationInputMatrix_t(populationInputMatrices,t):
    """
    Return the coastal matrix with population density at a given time, as interpolated
    from the set of populationInputMatrices
    """
    time_matrices = np.array([datetime(2000,1,1),datetime(2005,1,1),datetime(2010,1,1),
                               datetime(2015,1,1),datetime(2020,1,1)])
    if (t < time_matrices[0]) or (t>time_matrices[-1]):
        print('WARNING: specified time falls outside of range of population density datafiles: %s' % t)
    
    if t in time_matrices:
        arg = np.where(time_matrices == t)[0][0]
        populationMatrix = populationInputMatrices[arg,:,:]
    else:
        arg_t = time_matrices.searchsorted(t,side='left')
#        print(arg_t)

        dt_left = (t - time_matrices[arg_t-1]).total_seconds()
        dt_right = (time_matrices[arg_t] - t).total_seconds()
        
        a = dt_right/(dt_left+dt_right)
        b = dt_left/(dt_left+dt_right)
        
        # create linear interpolation of population densities
        populationMatrix = a*populationInputMatrices[arg_t-1,:,:] + b*populationInputMatrices[arg_t,:,:]
    

    return populationMatrix


def riverInputMatrix_t(riverInputMatrices,t,total_pop_t,scale_waste=True):
    """
    Return the matrix with river waste input at a given time, as given per month
    in the table from Lebreton 2017, corrected for population density
    
    total_pop_t: array with populations at given time; used to scale the river waste with, 
    assuming data from Lebreton is valid at 2010-01-01
    """
    
    t_month = t.month
    
    if scale_waste:
        total_pops = total_pop_t[0,:]
        time_pops = total_pop_t[1,:]
        
        if (t < time_pops[0]) or (t>time_pops[-1]):
            print('WARNING: specified time falls outside of range of population density datafiles: %s' % t)
        
        if t in time_pops:
            arg = np.where(time_pops == t)[0][0]
            totPop_t = total_pops[arg]
        else:
            arg_t = time_pops.searchsorted(t,side='left')
    #        print(arg_t)
    
            dt_left = (t - time_pops[arg_t-1]).total_seconds()
            dt_right = (time_pops[arg_t] - t).total_seconds()
            
            a = dt_right/(dt_left+dt_right)
            b = dt_left/(dt_left+dt_right)        
        
            totPop_t = a*total_pops[arg_t-1] + b*total_pops[arg_t]
        
        scaling = totPop_t/total_pops[2] #scale with total population in 2010   
    else:
        scaling = 1.0
        

    riverMatrix = scaling*riverInputMatrices[t_month-1,:,:]

    return riverMatrix


def calculateCoastalDistance(landMask,landBorderMask,fieldMesh_x,fieldMesh_y,do_plot=False):
    """
    Calculate distance to the closest coast in the ocean.
    Vectorized calculations, making it really fast
    """
    
    oceanCellIndices = np.where(landMask == 0)         #which indices in the grid correspond to ocean
    coastalCellIndices = np.where(landBorderMask == 1) #which correspond to coast (land)
    
    coastalDistance_mat = np.zeros([len(coastalCellIndices[0]),len(oceanCellIndices[0])]) #matrix in which distances are stored of given coastal loc. w.r.t. all ocean locations
    
    coastalDistance = np.zeros(fieldMesh_x.shape) #matrix corresponding to fieldMesh_x/y in which distances are stored
    
    for i1 in range(len(coastalCellIndices[0])): #go through all coastal elements
        lon_coast = fieldMesh_x[coastalCellIndices[0][i1],coastalCellIndices[1][i1]]    
        lat_coast = fieldMesh_y[coastalCellIndices[0][i1],coastalCellIndices[1][i1]]

        distMat_lon = (lon_coast - fieldMesh_x[oceanCellIndices[0],oceanCellIndices[1]])*1.11e2*0.788 #find distances coastal element w.r.t. ocean cells. 0.788 comes from lat=38deg (medit.)
        distMat_lat = (lat_coast - fieldMesh_y[oceanCellIndices[0],oceanCellIndices[1]])*1.11e2

        coastalDistance_mat[i1,:] = np.sqrt(np.power(distMat_lon, 2) + np.power(distMat_lat, 2))    
    

    minDist = np.min(coastalDistance_mat,axis=0) #from all coastal elements, get minimum distance to given ocean cell
    coastalDistance[oceanCellIndices[0],oceanCellIndices[1]] = minDist #store in matrix
    
    if do_plot:
        plt.figure()
        plt.contourf(fieldMesh_x,fieldMesh_y,coastalDistance,cmap=plt.cm.jet,levels=np.linspace(0,coastalDistance.max(),50))
        plt.colorbar()
        plt.title('Distance to closest coast [km]')
    
    return coastalDistance


def calculateLandCurrent(landMask,fieldMesh_x,fieldMesh_y,do_plot=False):
    """
    Calculate closest cell with water for all land cells, create vector field to
    closest water cell
    """
    oceanCellIndices = np.where(landMask == 0)         #which indices in the grid correspond to ocean
    landCellIndices = np.where(landMask == 1) #which correspond to land
        
    
    landVectorField_x = np.zeros(fieldMesh_x.shape)
    landVectorField_y = np.zeros(fieldMesh_y.shape)
    
    for i1 in range(len(landCellIndices[1])): #go through all land cells
        lon_coast = fieldMesh_x[landCellIndices[0][i1],landCellIndices[1][i1]]    
        lat_coast = fieldMesh_y[landCellIndices[0][i1],landCellIndices[1][i1]]

        distMat_lon = (lon_coast - fieldMesh_x[oceanCellIndices[0],oceanCellIndices[1]]) #find distances coastal element w.r.t. ocean cells. 
        distMat_lat = (lat_coast - fieldMesh_y[oceanCellIndices[0],oceanCellIndices[1]])

        distance_toOcean = np.sqrt(np.power(distMat_lon, 2) + np.power(distMat_lat, 2))    
        minDist = np.min(distance_toOcean)
        i_minDist = np.where(distance_toOcean == minDist)
        if len(i_minDist[0]) == 1:
            #easy case: vector to single point
            lon_ocean = fieldMesh_x[oceanCellIndices[0][i_minDist],oceanCellIndices[1][i_minDist]]
            lat_ocean = fieldMesh_y[oceanCellIndices[0][i_minDist],oceanCellIndices[1][i_minDist]]
            
            landVectorField_x[landCellIndices[0][i1],landCellIndices[1][i1]] = (lon_ocean - lon_coast) / np.sqrt((lon_ocean - lon_coast)**2 + (lat_ocean - lat_coast)**2)
            landVectorField_y[landCellIndices[0][i1],landCellIndices[1][i1]] = (lat_ocean - lat_coast) / np.sqrt((lon_ocean - lon_coast)**2 + (lat_ocean - lat_coast)**2)
        
        elif len(i_minDist[0]) > 1:
            #multiple ocean cells are the closest: take mean x and y vals
            lon_ocean = np.mean(fieldMesh_x[oceanCellIndices[0][i_minDist],oceanCellIndices[1][i_minDist]])
            lat_ocean = np.mean(fieldMesh_y[oceanCellIndices[0][i_minDist],oceanCellIndices[1][i_minDist]])
            
            landVectorField_x[landCellIndices[0][i1],landCellIndices[1][i1]] = (lon_ocean - lon_coast) / np.sqrt((lon_ocean - lon_coast)**2 + (lat_ocean - lat_coast)**2)
            landVectorField_y[landCellIndices[0][i1],landCellIndices[1][i1]] = (lat_ocean - lat_coast) / np.sqrt((lon_ocean - lon_coast)**2 + (lat_ocean - lat_coast)**2)            


    if do_plot:
        plt.figure()
        plt.quiver(landVectorField_x[::3,::3],landVectorField_y[::3,::3],scale=30)
    
    return landVectorField_x,landVectorField_y


def writeParticleInfo(fname,t,lat,lon,popDen,fishDen,t0=datetime(2001,1,1)):
    """
    Write release info particles for efficient processing
    """
    
    t_epoch = t-t0
    
    f = open(fname,'w')
    for i1 in range(len(t)):
        f.write('%i\t%3.10f\t%3.10f\t%3.10f\t%3.10f\n' % (t_epoch[i1].total_seconds(),lat[i1],lon[i1],popDen[i1],fishDen[i1]) )
    f.close()  


class PlasticParticle(JITParticle):
    age = Variable('age', dtype=np.float32, initial=0.)
    # beached : 0 sea, 1 beached, 2 after non-beach dyn, 3 after beach dyn, 4 please unbeach
    beached = Variable('beached',dtype=np.int32,initial=0.)
    
    coastalAge_1cell = Variable('coastalAge_1cell', dtype=np.float32, initial=0.)
    coastalAge_2cell = Variable('coastalAge_2cell', dtype=np.float32, initial=0.)
    coastalZoneAge   = Variable('coastalZoneAge', dtype=np.float32, initial=0.)
    
    dot_curr = Variable('dot_curr', dtype=np.float32, initial=0.)
    dot_stokes = Variable('dot_stokes', dtype=np.float32, initial=0.)
    dot_wind = Variable('dot_wind', dtype=np.float32, initial=0.)
    dot_wind_Ek = Variable('dot_wind_Ek', dtype=np.float32, initial=0.)
    
    mag_curr = Variable('mag_curr', dtype=np.float32, initial=0.)
    mag_stokes = Variable('mag_stokes', dtype=np.float32, initial=0.)
    mag_wind = Variable('mag_wind', dtype=np.float32, initial=0.)
    

def Ageing(particle, fieldset, time):
    dtt = particle.dt
    particle.age += dtt
    
def delete_180(particle, fieldset, time):
    if particle.age > (180*24*60*60):
        particle.delete()

def delete_300(particle, fieldset, time):
    if particle.age > (300*24*60*60):
        particle.delete()

def coastalAgeing(particle, fieldset, time):
    
    coastalDistance_p = fieldset.coastalDistance[time, particle.depth, particle.lat, particle.lon] 
    coastalZone_p = fieldset.coastalZone[time, particle.depth, particle.lat, particle.lon] 
    dtt = particle.dt
    
    if coastalDistance_p < 1.11e2/16:
        particle.coastalAge_1cell += dtt
    if coastalDistance_p < 2*(1.11e2/16):
        particle.coastalAge_2cell += dtt
    if (coastalZone_p > 0.001) and (coastalZone_p < 1):
        particle.coastalZoneAge += dtt


def coastalDynamics(particle, fieldset, time):
    coastalDistance_p = fieldset.coastalDistance[time, particle.depth, particle.lat, particle.lon] 
    dtt = particle.dt
        
    if coastalDistance_p < 2*(1.11e2/16):
        
        (u_curr, v_curr) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        (u_stokes, v_stokes) = fieldset.UV_Stokes[time, particle.depth, particle.lat, particle.lon]
        (u_wind, v_wind) = fieldset.UV_wind[time, particle.depth, particle.lat, particle.lon]
        
        cos_m45 = 0.7071067811865475
        sin_m45 = -0.7071067811865475
        
        u_wind_Ek = cos_m45*u_wind - sin_m45*v_wind
        v_wind_Ek = sin_m45*u_wind + cos_m45*v_wind
        
        lon_spacing = 0.0625#fieldset.U.lon[1] - fieldset.U.lon[0]
        lat_spacing = 0.0625#fieldset.U.lat[1] - fieldset.U.lat[0]
        
        land_top = fieldset.landMask[time, particle.depth, particle.lat+lat_spacing, particle.lon]
        
        land_bot = fieldset.landMask[time, particle.depth, particle.lat-lat_spacing, particle.lon]
        
        if particle.lon > -5.9375:
            land_left = fieldset.landMask[time, particle.depth, particle.lat, particle.lon-lon_spacing]
        else:
            land_left = 0.0 #TODO: this is specific for the medit., make more general
        
        land_right = fieldset.landMask[time, particle.depth, particle.lat, particle.lon+lon_spacing]
        
        dot_curr = 0.
        dot_stokes = 0.
        dot_wind = 0.
        dot_wind_Ek = 0.

        if land_top > 0.99:
            
            dot_curr_t = u_curr*0. + v_curr*1.
            dot_stokes_t = u_stokes*0. + v_stokes*1.
            dot_wind_t = u_wind*0. + v_wind*1.
            dot_wind_Ek_t = u_wind_Ek*0. + v_wind_Ek*1.
            
            if dot_curr_t > dot_curr:
                dot_curr = dot_curr_t 
            if dot_stokes_t > dot_stokes:
                dot_stokes = dot_stokes_t 
            if dot_wind_t > dot_wind:
                dot_wind = dot_wind_t 
            if dot_wind_Ek_t > dot_wind_Ek:
                dot_wind_Ek = dot_wind_Ek_t

        if land_bot > 0.99:
            
            dot_curr_t= u_curr*0. + v_curr*-1.
            dot_stokes_t = u_stokes*0. + v_stokes*-1.
            dot_wind_t = u_wind*0. + v_wind*-1.
            dot_wind_Ek_t = u_wind_Ek*0. + v_wind_Ek*-1.
            
            if dot_curr_t > dot_curr:
                dot_curr = dot_curr_t
            if dot_stokes_t > dot_stokes:
                dot_stokes = dot_stokes_t
            if dot_wind_t > dot_wind:
                dot_wind = dot_wind_t
            if dot_wind_Ek_t > dot_wind_Ek:
                dot_wind_Ek = dot_wind_Ek_t        
        
        if land_left > 0.99:
            
            dot_curr_t = u_curr*-1. + v_curr*0.
            dot_stokes_t = u_stokes*-1. + v_stokes*0.
            dot_wind_t = u_wind*-1. + v_wind*0.
            dot_wind_Ek_t = u_wind_Ek*-1. + v_wind_Ek*0.
            
            if dot_curr_t > dot_curr:
                dot_curr = dot_curr_t
            if dot_stokes_t > dot_stokes:
                dot_stokes = dot_stokes_t
            if dot_wind_t > dot_wind:
                dot_wind = dot_wind_t      
            if dot_wind_Ek_t > dot_wind_Ek:
                dot_wind_Ek = dot_wind_Ek_t        
        
        if land_right > 0.99:
            
            dot_curr_t = u_curr*1. + v_curr*0.
            dot_stokes_t = u_stokes*1. + v_stokes*0.
            dot_wind_t = u_wind*1. + v_wind*0.
            dot_wind_Ek_t = u_wind_Ek*1. + v_wind_Ek*0.
            
            if dot_curr_t > dot_curr:
                dot_curr = dot_curr_t
            if dot_stokes_t > dot_stokes:
                dot_stokes = dot_stokes_t
            if dot_wind_t > dot_wind:
                dot_wind = dot_wind_t
            if dot_wind_Ek_t > dot_wind_Ek:
                dot_wind_Ek = dot_wind_Ek_t

        particle.dot_curr += dtt*dot_curr
        particle.dot_stokes += dtt*dot_stokes
        particle.dot_wind += dtt*dot_wind
        particle.dot_wind_Ek += dtt*dot_wind_Ek
        
        particle.mag_curr += dtt*math.sqrt(u_curr**2+v_curr**2)
        particle.mag_wind += dtt*math.sqrt(u_wind**2+v_wind**2)
        particle.mag_stokes += dtt*math.sqrt(u_stokes**2+v_stokes**2)
        

def AdvectionRK4(particle, fieldset, time):
    if particle.beached == 0:
        (u1, v1) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
        
        (u2, v2) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat1, lon1]
        lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
        
        (u3, v3) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat2, lon2]
        lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
        
        (u4, v4) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3]
        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        
        particle.beached = 2


def StokesUV(particle, fieldset, time):
    if particle.beached == 0:
            dtt = particle.dt
            (u_uss, v_uss) = fieldset.UV_Stokes[time, particle.depth, particle.lat, particle.lon]

            particle.lon += u_uss * dtt
            particle.lat += v_uss * dtt
            particle.beached = 3


        
def WindUV(particle,fieldset,time):
    
    if particle.beached == 0:
        dtt = particle.dt
        (u_wnd, v_wnd) = fieldset.UV_wind[time, particle.depth, particle.lat, particle.lon]
        particle.lon += u_wnd * dtt
        particle.lat += v_wnd * dtt
        particle.beached = 3


def BrownianMotion2D(particle, fieldset, time):
    """Kernel for simple Brownian particle diffusion in zonal and meridional direction.
    Assumes that fieldset has fields Kh_zonal and Kh_meridional"""
    if particle.beached == 0:
        dtt = particle.dt
        r = 1/3.
        kh_meridional = fieldset.Kh_meridional[time, particle.depth, particle.lat, particle.lon]
        particle.lat += random.uniform(-1., 1.) * math.sqrt(2*math.fabs(dtt)*kh_meridional/r)
        kh_zonal = fieldset.Kh_zonal[time, particle.depth, particle.lat, particle.lon]
        particle.lon += random.uniform(-1., 1.) * math.sqrt(2*math.fabs(dtt)*kh_zonal/r)
        
        particle.beached = 3

def BeachTesting(particle, fieldset, time):
    if particle.beached == 2 or particle.beached == 3:
        (u, v) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        if u == 0 and v == 0:
            if particle.beached == 2:
                particle.beached = 4
            else:
                particle.beached = 1
        else:
            particle.beached = 0

def UnBeaching(particle, fieldset, time):
    if particle.beached == 4 or particle.beached == 1:
        dtt = particle.dt
        (u_land, v_land) = fieldset.UV_unbeach[time, particle.depth, particle.lat, particle.lon]
        particle.lon += u_land * dtt
        particle.lat += v_land * dtt
        particle.beached = 0

def deleteParticle(particle, fieldset, time):
    particle.delete()
       

 #%%   
if __name__ == "__main__":
    
    #%% Define variables
    K,fname,fraction_pop,fraction_fisheries,particles_per_day = getArguments(sys.argv[1:])
   
    fraction_rivers = 1 - (fraction_pop + fraction_fisheries)
    
    if fraction_rivers <= 0:
        sys.exit('fraction pop/fisheries should not add up to a value higher than 1')
    
    print('Running with K = %f, waste fraction pop = %f, waste fraction rivers = %f, waste fraction fisheries = %f, particles per day = %i' % (K,fraction_pop,fraction_rivers,fraction_fisheries,particles_per_day) )
        
    plottingDomain      = [23,38,30.5,38] #east mediterranean close-up
    day_start           = datetime(2006,2,1,0,0) #datetime(2001,1,1,0,0)    
    day_end             = datetime(2016,12,31,0,0) #datetime(2016,12,31,0,0)   
    dt_particle         = 'day'                            #dt for releasing particles, 'day' or 'week'
    dt_plot             = timedelta(days=14)
    projectFolder       = fname + '/'
    use_rivers          = True
    use_popDensity      = True
    use_fisheries       = True
    particlesReleaseInfo= 'particlesReleaseInfo.txt'    # data where and when particles are released for later analysis
    windage             = 0.01
    
    do_K                = True
    do_wind             = False
    do_stokes           = True
    do_unbeach          = True
   
    do_delete_180       = True
    do_delete_300       = False
    
    do_plot             = False
    

    n_plots         = int(np.floor((day_end-day_start)/dt_plot))
    if n_plots > 60:
        print('Creating %i plots and datafiles' % n_plots)

    if do_delete_180:
        print('Deleting particles at 180 days')
    elif do_delete_300:
        print('Deleting particles at 300 days')
    else:
        print('not deleting particles')
        
    #%% Set folders
   
    if os.environ['USER'] == 'kaandorp': # desktop
        homeDir = os.environ['HOME']
        outDir = os.path.join(os.getcwd(), projectFolder)
    

    if os.path.exists(outDir):
        print('Writing files to %s\n' % outDir)
    else:
        os.makedirs(outDir)
        print('Creating folder %s for output files\n' % outDir)
        
        

    
    #%% Set-up parcels fields and initialize particles/currents
    
    #----------------------- Ocean Currents -----------------------------------
    fileList = list(np.sort(glob.glob(os.path.join(homeDir , 'Data/CMEMS/Mediterranean/*'))))
    file_masks = fileList[0]
    mask_field = 'vozocrtx'

    
    filenames = {'U': fileList,
                 'V': fileList}
    variables = {'U': 'vozocrtx',
                 'V': 'vomecrty'}
    dimensions = {'lat': 'lat',
                  'lon': 'lon',
                  'time': 'time'}
    
    
    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions)
    lons = fieldset.U.lon
    lats = fieldset.U.lat
    fieldMesh_x,fieldMesh_y = np.meshgrid(lons,lats)
    
    
    #----------------------------- Wind ---------------------------------------
    fileList_wnd = list(np.sort(glob.glob(os.path.join(homeDir , 'Data/Wind/WW3/**/*'))))    
    
    filenames_wnd = {'U_wind': fileList_wnd,
                     'V_wind': fileList_wnd}
    variables_wnd = {'U_wind': 'uwnd',
                     'V_wind': 'vwnd'}
    dimensions_wnd = {'lat': 'latitude',
                      'lon': 'longitude',
                      'time': 'time'}    

    fieldset_wnd = FieldSet.from_netcdf(filenames_wnd, variables_wnd, dimensions_wnd)
    fieldset_wnd.U_wind.units = GeographicPolar()
    fieldset_wnd.V_wind.units = Geographic()
    fieldset_wnd.U_wind.set_scaling_factor(windage)
    fieldset_wnd.V_wind.set_scaling_factor(windage)

    fieldset.add_field(fieldset_wnd.U_wind)
    fieldset.add_field(fieldset_wnd.V_wind)
    
    vectorField_wnd = VectorField('UV_wind',fieldset.U_wind,fieldset.V_wind)
    fieldset.add_vector_field(vectorField_wnd)


    #----------------------------- Stokes -------------------------------------
    fileList_S_x = list(np.sort(glob.glob(os.path.join(homeDir , 'Data/CMEMS/Mediterranean_waves/VSDX/*'))))    
    fileList_S_y = list(np.sort(glob.glob(os.path.join(homeDir , 'Data/CMEMS/Mediterranean_waves/VSDY/*'))))    
    
    filenames_S = {'U_stokes': fileList_S_x, #Cannot be U for codegenerator!!
                   'V_stokes': fileList_S_y}
    variables_S = {'U_stokes': 'VSDX',
                   'V_stokes': 'VSDY'}
    dimensions_S = {'lat': 'lat',
                    'lon': 'lon',
                    'time': 'time'}    

    fieldset_Stokes = FieldSet.from_netcdf(filenames_S, variables_S, dimensions_S,mesh='spherical')
    fieldset_Stokes.U_stokes.units = GeographicPolar()
    fieldset_Stokes.V_stokes.units = Geographic()
    
    fieldset.add_field(fieldset_Stokes.U_stokes)
    fieldset.add_field(fieldset_Stokes.V_stokes)
    
    vectorField_Stokes = VectorField('UV_Stokes',fieldset.U_stokes,fieldset.V_stokes)
    fieldset.add_vector_field(vectorField_Stokes)
        

    #--------------- add land current, removing parcels from land -------------
    print('Calculating land/coast masks...')
    file_landMask = os.path.join('.',('datafile_landMask_%ix_%iy' % (len(lons),len(lats)) ) )
    file_coastMask = os.path.join('.',('datafile_coastMask_%ix_%iy' % (len(lons),len(lats)) ) )
    file_landBorderMask = os.path.join('.',('datafile_landBorderMask_%ix_%iy' % (len(lons),len(lats)) ) )
    file_coastalDistance = os.path.join('.',('datafile_coastalDistance_%ix_%iy' % (len(lons),len(lats)) ) )
    
    if os.path.exists(file_landMask):
        landMask = np.loadtxt(file_landMask)
    else:
        landMask = getLandMask(file_masks,mask_field) 
        np.savetxt(file_landMask,landMask)
        
    if os.path.exists(file_coastMask):
        coastMask = np.loadtxt(file_coastMask)
    else:
        coastMask = getCoastMask(landMask,lons,lats)
        np.savetxt(file_coastMask,coastMask)

    if os.path.exists(file_landBorderMask):
        landBorderMask = np.loadtxt(file_landBorderMask)
    else:
        landBorderMask = getLandBorder(landMask,lons,lats)
        np.savetxt(file_landBorderMask,landBorderMask)
        
    if os.path.exists(file_coastalDistance):
        coastalDistance = np.loadtxt(file_coastalDistance)
    else:
        coastalDistance = calculateCoastalDistance(landMask,landBorderMask,fieldMesh_x,fieldMesh_y)
        np.savetxt(file_coastalDistance,coastalDistance)
        
    #----------------------------- Unbeaching current -------------------------------------
    if do_unbeach:

        file_landcurrent_U = os.path.join('.',('datafile_landCurrentU_%ix_%iy' % (len(lons),len(lats)) ) )
        file_landcurrent_V = os.path.join('.',('datafile_landCurrentV_%ix_%iy' % (len(lons),len(lats)) ) )
        
        if os.path.exists(file_landcurrent_U):
            print('Using unbeaching, reading in landcurrent files...')
            landCurrent_U = np.loadtxt(file_landcurrent_U)
            landCurrent_V = np.loadtxt(file_landcurrent_V)
        else:
            landCurrent_U,landCurrent_V = calculateLandCurrent(landMask,fieldMesh_x,fieldMesh_y)
            np.savetxt(file_landcurrent_U,landCurrent_U)
            np.savetxt(file_landcurrent_V,landCurrent_V)
        
        U_land = Field('U_land',landCurrent_U,lon=lons,lat=lats,fieldtype='U',mesh='spherical')
        V_land = Field('V_land',landCurrent_V,lon=lons,lat=lats,fieldtype='V',mesh='spherical')
        
        fieldset.add_field(U_land)
        fieldset.add_field(V_land)
        
        vectorField_unbeach = VectorField('UV_unbeach',U_land,V_land)
        fieldset.add_vector_field(vectorField_unbeach)
    
    
    
    K_m = K*np.ones(fieldMesh_x.shape)
    K_z = K*np.ones(fieldMesh_x.shape)


    Kh_meridional = Field('Kh_meridional', K_m,lon=lons,lat=lats,mesh='spherical')
    Kh_zonal = Field('Kh_zonal', K_z,lon=lons,lat=lats,mesh='spherical')
    coastalDistance = Field('coastalDistance',coastalDistance,lon=lons,lat=lats,mesh='spherical')
    coastalZone = Field('coastalZone',coastMask,lon=lons,lat=lats,mesh='spherical')
    
    fieldset.add_field(Kh_meridional)
    fieldset.add_field(Kh_zonal)
    fieldset.add_field(coastalDistance)
    fieldset.add_field( Field('landMask',landMask,lon=lons,lat=lats,mesh='spherical') )
    fieldset.add_field(coastalZone)
    
    #%% Read and initialize river pollution data
    if use_rivers:
        file_riverInputMatrix = os.path.join('.',('datafile_riverInputMatrix_%ix_%iy' % (len(lons),len(lats)) ) )
        
        if os.path.exists(file_riverInputMatrix):
            # read file
            print('Reading in river input matrix files...')
            matrix_filefmt = np.loadtxt(file_riverInputMatrix )
            riverInputMatrices = np.reshape(matrix_filefmt,[12,int(matrix_filefmt.shape[0]/12),matrix_filefmt.shape[1]])
            
        else:
            rivers_shapefile    = os.path.join(homeDir , 'Data/PlasticData/PlasticRiverInputs_Lebreton/PlasticRiverInputs.shp')
            rivers_waste        = os.path.join(homeDir , 'Data/PlasticData/PlasticRiverInputs_Lebreton/PlasticRiverInputs.csv')
            riverInputMatrices  = initializeRiverMatrices(rivers_shapefile, rivers_waste, coastMask, fieldMesh_x, fieldMesh_y, selection_estim='mid')
            print('Saving river input matrix files: %s...' % file_riverInputMatrix)
            matrix_filefmt = riverInputMatrices.reshape([12*riverInputMatrices.shape[1],riverInputMatrices.shape[2]])
            np.savetxt(file_riverInputMatrix,matrix_filefmt)
    else:
        riverInputMatrices = np.zeros([12,fieldMesh_x.shape[0],fieldMesh_x.shape[1]])


    #%% Initialize pollution proportional to population density
    if use_popDensity:
        file_popInputMatrix = os.path.join('.',('datafile_popInputMatrix_50km_%ix_%iy' % (len(lons),len(lats)) ) )
        
        if os.path.exists(file_popInputMatrix):
            # read file
            print('Reading in population input matrix files...')
            matrix_filefmt = np.loadtxt(file_popInputMatrix )
            populationInputMatrices = np.reshape(matrix_filefmt,[5,int(matrix_filefmt.shape[0]/5),matrix_filefmt.shape[1]])
        else:
            pop_datafile            = os.path.join(homeDir , 'Data/populationDensity/gpw_v4_e_atotpopbt_dens_2pt5_min.nc')
            populationInputMatrices = initializePopulationMatrices_50km(pop_datafile,coastMask,fieldMesh_x,fieldMesh_y,filter_width=2,plotPopData=False)
            print('Saving population input matrix files: %s...' % file_popInputMatrix)
            matrix_filefmt = populationInputMatrices.reshape([5*populationInputMatrices.shape[1],populationInputMatrices.shape[2]])
            np.savetxt(file_popInputMatrix,matrix_filefmt)

    else:
        populationInputMatrices = np.zeros([5,fieldMesh_x.shape[0],fieldMesh_x.shape[1]])
        # years: 2000,2005,2010,2015,2020
    
    # array containing total coastal population at a given time    
    total_pop_t = np.array([np.sum(populationInputMatrices,axis=(1,2)),np.array([datetime(2000,1,1),datetime(2005,1,1),datetime(2010,1,1),
                               datetime(2015,1,1),datetime(2020,1,1)])])

    

    #%% Initialize plastic input from fisheries   
    if use_fisheries:
        file_fisheriesMatrix = os.path.join('.',('datafile_fisheries_%ix_%iy' % (len(lons),len(lats)) ) )
        
        if os.path.exists(file_fisheriesMatrix):
            print('Reading in fisheries input matrix files...')
            fisheriesInputMatrices = np.loadtxt(file_fisheriesMatrix )
            fisheriesInputMatrices[fisheriesInputMatrices<np.quantile(fisheriesInputMatrices[fisheriesInputMatrices>0],0.25)] = 0
#            fisheriesInputMatrices = np.reshape(matrix_filefmt,[12,int(matrix_filefmt.shape[0]/12),matrix_filefmt.shape[1]])
        else:
            print('error: fishery density matrix file not found: %s' % file_fisheriesMatrix)
            
    else:
        fisheriesInputMatrices = np.zeros([fieldMesh_x.shape[0],fieldMesh_x.shape[1]])   
        
    #%% initialize other stuff

    
 
    if dt_particle == 'day':
        n_days = (day_end-day_start).days
        releaseTimes = np.arange(0, n_days) * timedelta(hours=24)
    else:
        print('dt not implemented, releasing particles at once')
        releaseTimes = 0
    releaseTimes += day_start 
 
    
    # calculate waste for year of calibration (2010)
    # calibrate waste proportions, to have a yearly proportion according to scaling_popToRiver
    riverWaste_calib = 0
    popWaste_calib = 0
    for i1 in range(12):
        t_calibrate = datetime(2010,1+i1,1) # time at which proportion river/coastal waste is calibrated
        riverWaste_calib += np.sum(riverInputMatrix_t(riverInputMatrices,t_calibrate,total_pop_t))
        popWaste_calib += np.sum(populationInputMatrix_t(populationInputMatrices,t_calibrate))
    
    #input which doesn't vary monthly: fisheries
    fisheriesWaste_calib = 12*np.sum(fisheriesInputMatrices)
    
    # everything is scaled accouring to the river waste, since for this we have an idea how much
    # plastic enters the basin in tonnes, see Lebreton (2017)
    scaling_popToRiver = fraction_pop / fraction_rivers
    scaling_fisheriesToRiver = fraction_fisheries / fraction_rivers
    
    scaling_pop = scaling_popToRiver*(riverWaste_calib/popWaste_calib) #scaling for popWaste: to tonnes/month
    scaling_fisheries = scaling_fisheriesToRiver*(riverWaste_calib/fisheriesWaste_calib)
    
    riverWaste_calib /= 12 
    popWaste_calib = scaling_pop*(popWaste_calib/12)
    fisheriesWaste_calib = scaling_fisheries*(fisheriesWaste_calib/12) # not necessary, but for consistency..
    
    wastePerMonth_calib = riverWaste_calib + popWaste_calib + fisheriesWaste_calib #average monthly waste used for calibration
    particle_weight = (12*wastePerMonth_calib)/(particles_per_day*365)
    print('Each particle stands for approx. %f tonnes of emitted plastic' % particle_weight)



    # -------- create figures for validation ----------------------------------    
    times_plot = np.array([])
    waste_river = np.array([])
    waste_pop = np.array([])
    waste_fisheries = np.array([])
    year_0 = day_start.year
    for i1 in range((day_end.year - day_start.year)+1):
        for i2 in range(12):
            t_plot = datetime(year_0+i1,i2+1,1)
            times_plot = np.append(times_plot,t_plot)
            waste_river = np.append(waste_river, np.sum(riverInputMatrix_t(riverInputMatrices,t_plot,total_pop_t)))
            waste_pop = np.append(waste_pop, np.sum(populationInputMatrix_t(populationInputMatrices,t_plot)))
            waste_fisheries = np.append(waste_fisheries,np.sum(fisheriesInputMatrices) )
    
    fig, ax1 = plt.subplots(3,figsize=(20,16))
    ax1[0].plot_date(times_plot,waste_river, 'b-')
    ax1[0].set_ylabel('Total waste emitted by\n rivers [tonnes/month]', color='b')
    ax1[0].tick_params('y', colors='b')
    
    ax2 = ax1[0].twinx()
    ax2.plot_date(times_plot,waste_pop, 'r.')
    ax2.set_ylabel('Total population density\n at coast', color='r')
    ax2.tick_params('y', colors='r')
    
    ax3 = ax1[1].twinx()
    ax3.plot_date(times_plot,waste_fisheries, 'r.')
    ax3.set_ylabel('Fishery intensity', color='r')
    ax3.tick_params('y', colors='r')    
    
    ax1[2].plot_date(times_plot,waste_river + scaling_pop*waste_pop + scaling_fisheries*waste_fisheries, 'k-')
    ax1[2].set_ylabel('Total modelled waste emitted,\n [tonnes/month]')
    ax1[2].set_xlabel('date')

#    ax3 = ax1[1].twinx()
    

    #%% Create arrays with particles for parcels

    particles_lat   = np.array([])
    particles_lon   = np.array([])
    particles_t     = np.array([])

    # arrays containing information on the release location
    particles_origPopDen = np.array([])
    particles_origFishDen = np.array([])

    # every time step, $total_particles are released, each particle has a probability
    # adding up to one to end up on one of the coast cells
    
    print('randomly drawing particles with probability proportional to the estimated input waste at coast, this can take a while...')
    n_particles_released = []
    t_start = time.time()
    for i1 in range(len(releaseTimes)):

        releaseTime = releaseTimes[i1]
        t_year = releaseTime.year
        t_month = releaseTime.month
        t_day = releaseTime.day
        
        _,n_days = calendar.monthrange(t_year,t_month)
        
        correction_ndays = (365/12)/n_days # correct for longer/shorter months compared to the average amount of days/month in 2010

        
        # get monthly waste [tonnes] at given time
        riverInput = correction_ndays*riverInputMatrix_t(riverInputMatrices,releaseTime,total_pop_t)
        popInput = scaling_pop*populationInputMatrix_t(populationInputMatrices,releaseTime)
        fisheriesInput = scaling_fisheries*fisheriesInputMatrices
        totInput = riverInput+popInput+fisheriesInput
        
        n_particles = particles_per_day*(np.sum(totInput)/wastePerMonth_calib) #*correction_ndays
        
        # convert n_particles, which is a float, to integer amount of particles, where the fraction
        # is converted to 0 or 1 depending on its magnitude        
        total_particles_base = np.floor(n_particles)
        chance_0 = 1-(n_particles-total_particles_base) # chances of adding 0/1
        chance_1 = n_particles-total_particles_base
    
        total_particles = int(total_particles_base + np.random.choice([0,1],size=1,p=[chance_0,chance_1])[0])
        
        
        # array with amount of particles released for later analysis
        n_particles_released.append(total_particles)
            
        # probability matrix of where the particle is released (sum=1)        
        totInput_prob = totInput / np.sum(totInput)
        probArray = totInput_prob[totInput_prob>0]
    
        #total_particles ~ total waste generated at certain time
        for i2 in range(total_particles):
    
            drawParticle = np.random.rand(len(probArray)) < probArray # adapted from  https://stackoverflow.com/questions/40474436/how-to-apply-numpy-random-choice-to-a-matrix-of-probability-values-vectorized-s
            lats_draw = fieldMesh_y[totInput_prob>0][drawParticle]
            lons_draw = fieldMesh_x[totInput_prob>0][drawParticle]
    
            particles_lat = np.append(particles_lat,lats_draw)
            particles_lon = np.append(particles_lon,lons_draw)
            particles_t = np.append(particles_t,np.array([releaseTime for i in range(len(lons_draw))]))   

        
        if (i1 % (int(len(releaseTimes)/10))) == 0:
            t_curr = time.time()
            print('%i / %i, %3.3f minutes' % (i1,len(releaseTimes), (t_curr-t_start)/60))
        
    #save particle origin popdensity and tourismdensity
    for lat,lon,t in zip(particles_lat,particles_lon,particles_t):
        popInput = populationInputMatrix_t(populationInputMatrices,t)

        index_lat = np.where(lats == lat)[0][0]
        index_lon = np.where(lons == lon)[0][0]
        particles_origPopDen = np.append(particles_origPopDen,popInput[index_lat,index_lon])
        particles_origFishDen = np.append(particles_origFishDen,fisheriesInputMatrices[index_lat,index_lon])


    ax3 = ax1[2].twinx()
    ax3.plot_date(releaseTimes,n_particles_released, 'g.')
    ax3.set_ylabel('n particles released', color='g')
    ax3.tick_params('y', colors='g')
    plt.savefig(os.path.join(outDir , 'ModelledWasteOutput.png') )
    
    writeParticleInfo(os.path.join(outDir,particlesReleaseInfo),particles_t,particles_lat,particles_lon,particles_origPopDen,particles_origFishDen,t0=day_start)
    print('Particle release info written to %s' % os.path.join(outDir,particlesReleaseInfo))
    
    

    pset = ParticleSet(fieldset, PlasticParticle, lon=particles_lon, lat=particles_lat, 
                       time=particles_t, repeatdt=None)


    kernel = pset.Kernel(AdvectionRK4) 
    if do_unbeach:
        kernel += pset.Kernel(BeachTesting) + pset.Kernel(UnBeaching)
    if do_K:
        kernel+= pset.Kernel(BrownianMotion2D) + pset.Kernel(BeachTesting) + pset.Kernel(UnBeaching)
    if do_wind:
        kernel += pset.Kernel(WindUV) + pset.Kernel(BeachTesting) + pset.Kernel(UnBeaching)
    if do_stokes:
        kernel += pset.Kernel(StokesUV) + pset.Kernel(BeachTesting) + pset.Kernel(UnBeaching)
    kernel += pset.Kernel(Ageing) + pset.Kernel(coastalAgeing) + pset.Kernel(coastalDynamics)
      
    if do_delete_180:
        kernel += pset.Kernel(delete_180)
    if do_delete_300:
        kernel += pset.Kernel(delete_300)
        
    #%% Run parcels
    for cnt in range(n_plots):
        plt.close('all')
        pset.execute(kernel,
                     runtime=dt_plot,  # runtime controls the interval of the plots
                     dt=timedelta(minutes=20),
                     recovery={ErrorCode.ErrorOutOfBounds: deleteParticle},
                     output_file=pset.ParticleFile(name=os.path.join(outDir +'K%3.1f_data_%3.3i.nc'%(K,cnt)), outputdt=timedelta(hours=24)) )
        
        if do_plot:
            currentTime = pd.to_datetime(pset.time_origin.time_origin) + timedelta(seconds=pset[0].time)
            plotParticles(pset,fieldMesh_x,fieldMesh_y,plottingDomain,landMask,pset[0].time)    
            plt.title(str(currentTime))
            plt.savefig(os.path.join(outDir + 'K%3.1f_particles_%3.3i.png' % (K,cnt)))
            plt.pause(0.01)

   
   
   