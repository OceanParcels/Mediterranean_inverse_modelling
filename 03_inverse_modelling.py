#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 15:30:46 2018
Most important script: optimize the parameters to match the observations as well as possible
@author: kaandorp
"""

from datetime import datetime, timedelta
import numpy as np
import glob
from scipy.stats import norm, linregress, pearsonr, spearmanr
from scipy.interpolate import griddata
import pandas as pd
import os
from weighted_KDE import gaussian_kde
import dask.dataframe as dd
from dask.delayed import delayed
import particleKDEWeights
import cartopy.crs as ccrs
import traceback
import pickle
import matplotlib.colors as mcolors

class linregress_2:
    """
    Linear regression without intercept
    """
    def __init__(self,x,y):
        a,_,_,_ = np.linalg.lstsq(x[:,np.newaxis], y, rcond=None)
        self.slope = a[0]
        self.intercept = 0


class linregress_log:
    """
    Regression of a powerlaw (i.e. linear on a log-scale)
    """
    def __init__(self,x,y):

        log10_x = np.log10(x)
        log10_y = np.log10(y)
        
        self.slope = 10**( (sum(log10_y) - sum(log10_x) ) / len(log10_x) )
        self.intercept = 0


def load(filename):
    pd_df = pd.read_csv(filename,dtype=dtype_pd,usecols=usecols_pd)
    return pd_df


def initMeasured():
    """
    Read in measurement data
    """
    measured = []
    lon_meas = []
    lat_meas = []
    study_ID = []
    kukulka_meas = []
    sizeRange_corr = []
    datetime_meas = np.array([])
    unit_meas = []
    basin_data = np.loadtxt('datafile_basins')
    mask_basin = (basin_data > 0)
    lons = np.linspace(-6,36.25,677)
    lats = np.linspace(30.1875,45.9375,253)
    fieldMesh_x,fieldMesh_y = np.meshgrid(lons,lats)
    
    filenames_meas = np.sort(glob.glob(os.path.join(csvFiles_measurements,'*measurement*.csv')))
    
    for i1 in range(len(filenames_meas)):
        filename_base = os.path.basename(filenames_meas[i1]) 
        year_meas = int(filename_base[0:4])
        month_meas = int(filename_base[4:6])
        day_meas = int(filename_base[6:8])
        datetime_meas = np.append(datetime_meas,datetime(year_meas,month_meas,day_meas) )
        data_meas = pd.read_csv(filenames_meas[i1])
        measured.append(data_meas.loc[0,'measured density'])
        lon_meas.append(data_meas.loc[0,'lon'])
        lat_meas.append(data_meas.loc[0,'lat'])
        study_ID.append(data_meas.loc[0,'# study_ID'])
        kukulka_meas.append(data_meas.loc[0,'Kukulka factor'])
        sizeRange_corr.append(data_meas.loc[0,'size range factor'])
        
        if 'weight' in filenames_meas[i1]:
            unit_meas.append('g/km2')
        elif 'num' in filenames_meas[i1]:
            unit_meas.append('#/km2')
        else:
            raise RuntimeError('unknown unit')
                              
    basin = griddata((fieldMesh_x[mask_basin],fieldMesh_y[mask_basin]),basin_data[mask_basin],(np.array(lon_meas), np.array(lat_meas)),method='nearest')
    return np.array(measured), np.array(lon_meas), np.array(lat_meas), np.array(study_ID), basin, np.array(kukulka_meas), datetime_meas, np.array(unit_meas), np.array(sizeRange_corr)


#%% 

class forward_model(object):
    """
    forward model, calculating modelled output corresponding to given measurements, for a set of parameters x
    parameters:
        
    x:
        model parameters, np.array
    x_normalized:
        bool, if True, x runs from 0 to 1, if not, the true value is used
    xaxis:
        'modelled' or 'measured': parameter defining whether model or measurements are put on x axis (modelled is standard)
    do_Kukulka:
        bool, use Kukulka correction
    density_type:
        set to 'KDE' (kernel density estimate)
    fn_type:
        regression function to use ('log_slopeonly' should be used, other options are 'linear' or 'linear_slopeonly')
    correlation_type:
        'spearman' or 'pearson'
    unitless:
        bool, if False use true particle mass, if true scaling for mass is done similar to abundance
    study_ID_use:
        np.array defining which measurement campaigns are used
    """ 
    
    
    def __init__(self,x,x_normalized=True,xaxis='modelled',do_Kukulka=True,density_type='KDE',
                 fn_type='log_slopeonly',correlation_type='spearman',
                 unitless=False,draw_errorbars=True,print_stats=False,study_ID_use=np.arange(1,21) ):
        self.x = x
        self.x_normalized = x_normalized
        self.xaxis = xaxis
        self.fn_type = fn_type
        self.do_Kukulka = do_Kukulka
        self.correlation_type = correlation_type
        self.unitless = unitless
        self.draw_errorbars = draw_errorbars
        self.density_type = density_type
        self.print_stats = print_stats
        self.study_ID_use = study_ID_use
        
    def calc_modelled_measured(self,x_total,units_):    
        """
        calculate modelled values for a set of measurements, given a certain unit
        """
        t_beach_tau,t_sink_tau,t_sink_rate,sink_init,ratio_popToRiver,ratio_fishToRiver,kde_kernel_bw,factor_river = x_total
        
        unitless_ = self.unitless
        if units_ == 'g/km2': #'g' in units_:
            if unitless_:
                use_kgs = False
            else:
                use_kgs = True
            
        elif units_ == '#/km2':
            use_kgs = False
            unitless_ = True
        else:
            print('error: units unknown')
            
        xaxis = self.xaxis #modelled measured
        
        modelled = []
    
        
        for i1 in range(len(filenames)):
            filename_base = os.path.basename(filenames[i1])
            year_ = int(filename_base[0:4])
            month_ = int(filename_base[4:6])
            day_ = int(filename_base[6:8])
        
            indices_measurements = np.where((datetime_meas == datetime(year_,month_,day_)) & (unit_meas == units_))[0]
        
            df_i = df.iloc[indices_0[i1]:indices_0[i1+1],:]
    
            
            lon_p = df_i.loc[:,'lon'].values
            lat_p = df_i.loc[:,'lat'].values
            t_age = df_i.loc[:,'age'].values
            t_coastal = df_i.loc[:,'coastalAge_1cell'].values

            source = init_source[df_i.loc[:,'# ID']]
            init_month_id_ = init_month_id[df_i.loc[:,'# ID']]
            init_lat_id_ = init_lat_id[df_i.loc[:,'# ID']]
            init_lon_id_ = init_lon_id[df_i.loc[:,'# ID']]                                        
                       
            if np.where(np.isnan(t_coastal))[0]:
                print('nan found')
                print(year_,month_,day_)
                break

            # calculate the abundance and mass inside the virtual particles
            weights = particleKDEWeights.calculateKDEWeights(t_age,t_coastal,source,t_beach_tau,t_sink_tau,t_sink_rate,sink_init,
                                                             ratio_popToRiver,ratio_fishToRiver,factor_river,init_month_id_,init_lon_id_,init_lat_id_,riverInput_low,riverInput_mid,riverInput_high,
                                                             n_pop=n_pop,n_river=n_river,n_fish=n_fish,kg_particle_river_mid=kg_particle_river_mid,use_kgs=use_kgs)

            if self.print_stats:
                if units_ == 'g/km2' and i1 == (len(filenames)-1):                
                    log_weights = np.log10(weights)
                    sigma_meas = np.sqrt(0.2201)
                    log_weights_perturbed = log_weights +  np.random.normal(0,sigma_meas,len(log_weights))
                    weights_perturbed = 10**log_weights_perturbed
                    
                    print('total mass norm.: %f' % weights.sum())
                    print('total mass pert.: %f' % weights_perturbed.sum())
                

            if self.density_type == 'KDE':
                density_kernel = gaussian_kde(([lon_p,lat_p]),bw_method=kde_kernel_bw,weights=weights)       
            else:
                raise RuntimeError('density type not defined')

            # loop through the available measurements for the given day & unit
            for index_meas in indices_measurements:
                x_measurement = lon_meas[index_meas]
                y_measurement = lat_meas[index_meas]
            
                # give the kernel density estimate. When particles are given in kg, this results in kg/km2 as output 
                kde_meas = density_kernel([x_measurement,y_measurement]) * (weights.sum()) # times sum of the weights since the KDE is normalized to 1
               
                if self.density_type == 'KDE':
                    modelled.append(list(kde_meas)[0])
                else:
                    modelled.append(kde_meas)
        
        modelled = np.array(modelled)
     

        mask_units_use = (unit_meas == units_)
        if self.study_ID_use is None:
            mask_model = (measured[mask_units_use] != 0) # remove modelled values where the measurement was 0
            mask_meas = ((unit_meas == units_) & (measured !=0)) # keep measured values for the appropriate unit, remove 0 measurements
        else:
            mask_study_ID_use = np.isin(study_ID,self.study_ID_use)
            mask_model = (measured[mask_units_use] != 0) & (mask_study_ID_use[mask_units_use])
            mask_meas = (unit_meas == units_) & (measured !=0) & (mask_study_ID_use)
        
        
        if xaxis == 'measured':
            y_ = modelled[mask_model]
            if self.do_Kukulka:
                x_ = measured[mask_meas]*kukulka_meas[mask_meas]*sizeRangeCorr_meas[mask_meas]
            else:
                x_ = measured[mask_meas]

                
        elif xaxis == 'modelled':
            x_ = modelled[mask_model]
            if self.do_Kukulka:
                y_ = measured[mask_meas]*kukulka_meas[mask_meas]*sizeRangeCorr_meas[mask_meas]
            else:
                y_ = measured[mask_meas]


        if np.isnan(y_.min()):
            print('NaN in measurements!!!!')
            print(y_)
            print(xaxis)
            
        
        if self.fn_type == 'linear':
            self.fit_fn = lambda x,y: linregress(x,y)
        elif self.fn_type == 'linear_slopeonly':
            self.fit_fn = lambda x,y: linregress_2(x,y)
        elif self.fn_type == 'log_slopeonly':
            self.fit_fn = lambda x,y: linregress_log(x,y)
        else:
            print('unknown fitting function')
            
        linregres = self.fit_fn(x_,y_)


        
        if unitless_:
            x_scaled = linregres.intercept + linregres.slope*x_
        else:
            x_scaled = x_
                
        return x_scaled, y_, mask_meas


    def to_unitfull(self,x_unitless_,x_lower_,x_upper_,log_scale_):
        vals_unitfull_ = np.zeros(x_lower_.shape)
        for i1,bool_log in enumerate(log_scale_):
            if bool_log:
                vals_unitfull_[i1] = 10**(x_lower_[i1]+x_unitless_[i1]*(x_upper_[i1]-x_lower_[i1]))
            else:
                vals_unitfull_[i1] = x_lower_[i1]+x_unitless_[i1]*(x_upper_[i1]-x_lower_[i1])
        return vals_unitfull_

    def to_unitless(self,x_unitfull_,x_lower_,x_upper_,log_scale_):
        vals_unitless_ = np.zeros(x_lower_.shape)
        for i1,bool_log in enumerate(log_scale_):
            if bool_log:
                exponent = np.log10(x_unitfull_[i1])
                vals_unitless_[i1] = (exponent - x_lower_[i1]) / (x_upper_[i1]-x_lower_[i1])
            else:
                vals_unitless_[i1] = (x_unitfull_[i1] - x_lower_[i1]) / (x_upper_[i1]-x_lower_[i1])
        return vals_unitless_
        
    def get_parameter_vals(self):
        
        if self.x_normalized:
            # need to convert to unitfull parameters
            x_total = np.zeros(x_lower.shape)
            
            x_total[x_use] = self.to_unitfull(self.x[x_use],x_lower[x_use],x_upper[x_use],log_scale[x_use])
            x_total[~x_use] = x_standard[~x_use]
            
            x_unitless = np.zeros(x_lower.shape)
            x_unitless[x_use] = self.x[x_use]
            x_unitless[~x_use] = self.to_unitless(x_standard[~x_use],x_lower[~x_use],x_upper[~x_use],log_scale[~x_use])
            
        else:
            # convert to unitless afterwards for calculating mismatch w.r.t. a-priori parameter choice
            x_total = np.copy(x_standard)
            x_total[x_use] = self.x[x_use]
            
            x_unitless = self.to_unitless(x_total,x_lower,x_upper,log_scale)
            
        self.x_unitless = x_unitless
        
        return x_total


    def calc_fit_numweight(self,units_):
        
        x_total = self.get_parameter_vals()
        
        if self.correlation_type == 'spearman':
            self.corr_fn = lambda x,y: spearmanr(x,y)[0]
        elif self.correlation_type == 'pearson':
            self.corr_fn = lambda x,y: pearsonr(x,y)[0]
        
        success_ = False
        
        try:
            if units_ == 'g/km2': 
                i_u = 1
            elif units_ == '#/km2':
                i_u = 0
                
            x_,y_,mask_meas_use = self.calc_modelled_measured(x_total,units_)
            correlation = self.corr_fn(x_,y_)
            
            log10_modelled = np.log10(x_[(x_>0)&(y_>0)])
            log10_measured = np.log10(y_[(x_>0)&(y_>0)])

            
            if self.print_stats:
                # ------------------- print statistics ------------------------
                print('model; mean, std: %f, %f' % (log10_modelled.mean(),log10_modelled.std()))
               
                #calculate perturbed state, calculate mean/std
                log10_modelled_perturbed_m = np.array([(log10_modelled + np.random.normal(0,std_array[i_u],len(log10_modelled))).mean() for i in range(20)])
                log10_modelled_perturbed_s = np.array([(log10_modelled + np.random.normal(0,std_array[i_u],len(log10_modelled))).std() for i in range(20)])
                print('model perturbed, mean, std: %f, %f' % (log10_modelled_perturbed_m.mean(),log10_modelled_perturbed_s.mean()))
                print('meas.; mean, std: %f, %f' % (log10_measured.mean(),log10_measured.std()))


            RMSE = np.sqrt(np.mean((log10_modelled-log10_measured)**2))
            
            if np.isnan(RMSE):
                success_ = False
            else:
                success_ = True
                
        except Exception: 
            print('Exception occurred while calculating the modelled particle concentrations:')
            traceback.print_exc()
            x_ = np.array([])
            y_ = np.array([])
            log10_modelled = np.array([])
            log10_measured = np.array([])
            
            mask_meas_use = np.array([])
            RMSE = np.nan
            correlation = np.nan
            success_ = False
    
        if WRITE:
            f = open(os.path.join(outDir,optimFile), 'a')
            f.write(np.array_repr(x_total).replace('\n', '') + ('\t%f\t%f\n' % (RMSE,correlation)) )
            f.close()    
        

        self.x_ = x_
        self.y_ = y_
        self.mask_meas_use = mask_meas_use
        self.correlation = correlation
        self.units = units_

        return log10_modelled, log10_measured, RMSE, correlation, success_


    def calculate_matrices(self):
        """
        calculate necessary matrices for evaluating the fitness
        """
        assert(type(units) == list)
        all_success = True
        
        C_d = np.array([])
        g_m = np.array([])
        d   = np.array([])
        mask_measurements = []
        measurements = []
        modelled = []
        
        # surface measurements
        for i1,unit in enumerate(units):
            log10_modelled, log10_measured, RMSE, correlation, success_ = self.calc_fit_numweight(unit)
#            print(RMSE,success_)
            if not success_:
                all_success = False
    
            g_m = np.append(g_m, log10_modelled)
            d = np.append(d, log10_measured)
            C_d = np.append(C_d, (std_array[i1]**2)*np.ones(len(log10_modelled)) ) #variance: square the std
            mask_measurements.append(self.mask_meas_use)
            measurements.append(log10_measured)
            modelled.append(log10_modelled)


#        print(measurements)
#        print(modelled)
        C_d = np.diag(C_d)
        
        param_deviation = self.x_unitless - 0.5*np.ones(len(self.x_unitless))
        C_m = (1./6.)**2 * np.eye(len(self.x_unitless)) #1 std is set at 0.25 (2std = 1)
        
        return(g_m,d,C_d,C_m,param_deviation,all_success,mask_measurements,measurements,modelled)


    def fitness(self):
        """
        Fitness function
        """
        g_m,d,C_d,C_m,param_deviation,all_success = self.calculate_matrices()
        
        mismatch = g_m - d
 
        if all_success:
            fitness = np.dot(mismatch.T, np.dot(np.linalg.inv(C_d), mismatch)) + np.dot(param_deviation.T, np.dot(np.linalg.inv(C_m), param_deviation))
        else:
            fitness = 1e10
           
        if WRITE:
            f = open(os.path.join(outDir,optimFile), 'a')
            f.write('%f\n' % fitness)
            f.close()    
        
        return fitness


    def plot_correlation(self):
        """
        Wrapper for correlation plot, can accept list of units
        """        
        if type(units) == list:
            for unit in units:
                self.plot_correlation_unit(unit)
        elif type(units) == str:
            self.plot_correlation_unit(units)
        

    def plot_correlation2(self):
        """
        Wrapper for correlation plot, put everything in one neat plot
        """        
        if type(units) == list:
            fig, axes = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1,1]},figsize=(8.3,6))
            
            for unit in units:
                self.plot_correlation_unit2(unit,fig,axes)
                self.plot_mismatch_map(unit)
        
        elif type(units) == str:
            print('errorrrr')
        
        
    def plot_correlation_unit(self,units_):
        """
        Create correlation plot model vs measurements
        """
        self.calc_fit_numweight(units_)
        
        if units_ == '#/km2':
            units_string = r'#/km$^2$'
        elif units_ == 'g/km2':
            units_string = r'g/km$^2$'
        else:
            raise RuntimeError('unknown units')
    
        
        df_keys = np.array(['Baini et al., 2018','Collignon et al., 2012','Fossi et al., 2012','Galgani et al. 2011','Galgani et al. 2012','van der Hal et al., 2017','Panti et al., 2015',
                     u'Güven et al., 2017',u'Gajšt et al. 2016','Suaria et al., 2016','Cozar et al., 2015','Collignon et al., 2014','Suaria et al., 2014',
                     u'Ruiz-Orejón et al., 2016',u'Ruiz-Orejón et al., 2018',u'Gündoğdu, 2017',u'Gündoğdu, 2018',u'Gündoğdu, 2017','Zeri et al., 2018','Pedrotti et al., 2016','de Haan et al., 2019'])
    
                
        not_seen = np.ones(df_keys.shape,dtype=bool)    
        cm = plt.cm.tab20

        fig1, axes1 = plt.subplots(figsize=(7,6))
        for i1 in range(len(self.y_)):
            studyID_plt = int(study_ID[self.mask_meas_use][i1]) 
            if not_seen[studyID_plt]:
                axes1.plot(self.x_[i1],self.y_[i1],'o', color=cm(studyID_plt), label=df_keys[studyID_plt] )
                not_seen[studyID_plt] = False
            else:
                axes1.plot(self.x_[i1],self.y_[i1],'o', color=cm(studyID_plt) )
        
        if self.xaxis == 'modelled':
            axes1.set_ylabel(r'Measured plastic density [%s]' % units_string)
            axes1.set_xlabel('Modelled plastic density [%s]' % units_string) 
            plt.title('Measured vs. modelled plastic densities, R: %f' % self.correlation)
        else:
            axes1.set_xlabel(r'Measured plastic density [%s]' % units_string)
            axes1.set_ylabel('Modelled plastic density [%s]' % units_string) 
            plt.title('Modelled vs. measured plastic densities, R: %f' % self.correlation)            
        
        try:
            minval = 0.9*min(min(self.x_[~(self.x_ == 0)]),min(self.y_[~(self.y_ == 0)]))
        except:
            minval = 1
            
        try:
            maxval = 1.11*max(max(self.x_),max(self.y_))
        except:
            maxval = 1e7
            
        axes1.plot([minval,maxval],[minval,maxval],'k-',label='1:1')
        
    
        if self.draw_errorbars:
            n_std = 2
            if units_ == '#/km2':
                log10_std_y = std_array[0]

            elif units_ == 'g/km2':
                log10_std_y = std_array[1]

            dy = n_std*log10_std_y
            
            y1u = minval+(minval*10**dy-minval)
            y2u = maxval+(maxval*10**dy-maxval)
            
            y1l = minval-(minval-10**(-dy)*minval)
            y2l = maxval-(maxval-10**(-dy)*maxval)
            
            axes1.plot([minval,maxval],[y1u,y2u],'r-',label='2x std from variogram')
            axes1.plot([minval,maxval],[y1l,y2l],'r-')
            
            # calculate % within errorbars
            mat = np.array([[minval,1],[maxval,1]])
            a_u,b_u = np.dot( np.linalg.inv(mat),np.array([[y1u],[y2u]]) )
            a_l,b_l = np.dot( np.linalg.inv(mat),np.array([[y1l],[y2l]]) )
            
            data_upper = a_u[0]*self.x_ + b_u[0]
            data_lower = a_l[0]*self.x_ + b_l[0]

            data_within = (self.y_ < data_upper) & (self.y_ > data_lower)
            percentage_within = 100*(np.sum(data_within)/len(data_within))

            plt.title('Measured vs. modelled plastic densities, R: %f \n percentage within %i$\sigma$: %3.1f' % (self.correlation,n_std,percentage_within) )
            
        axes1.set_yscale('log')
        axes1.set_xscale('log')
        plt.legend()            

    def plot_mismatch_map(self,units_):
       
        log10_modelled, log10_measured, RMSE, correlation, success_ = self.calc_fit_numweight(units_)

        diff = log10_modelled - log10_measured
        max_abs_diff = np.max(np.abs(diff))
        plt.figure(figsize=(20,6.5))     
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines(resolution='10m',color='grey')
        plt_scatter = ax.scatter(lon_meas[self.mask_meas_use],lat_meas[self.mask_meas_use],c=diff,cmap=plt.cm.RdBu_r,vmin=-max_abs_diff,vmax=max_abs_diff)
        ax.set_title(units_)
        plt.colorbar(plt_scatter)
    
    
    def plot_correlation_unit2(self,units_,fig,axes):
        """
        Create correlation plot model vs measurements
        """
        self.calc_fit_numweight(units_)
        
        if units_ == '#/km2':
            units_string = r'#/km$^2$'
            index_ax = 0
        elif units_ == 'g/km2':
            units_string = r'g/km$^2$'
            index_ax = 1
        else:
            raise RuntimeError('unknown units')
    
        
        df_keys = np.array(['Baini et al., 2018','Collignon et al., 2012','Fossi et al., 2012','Galgani et al. 2011','Galgani et al. 2012','van der Hal et al., 2017','Panti et al., 2015',
                     u'Güven et al., 2017',u'Gajšt et al. 2016','Suaria et al., 2016','Cozar et al., 2015','Collignon et al., 2014','Suaria et al., 2014',
                     u'Ruiz-Orejón et al., 2016',u'Ruiz-Orejón et al., 2018',u'Gündoğdu, 2017',u'Gündoğdu, 2018',u'Gündoğdu, 2017','Zeri et al., 2018','Pedrotti et al., 2016','de Haan et al., 2019'])
    
                
        not_seen = np.ones(df_keys.shape,dtype=bool)           
        
        colors1 = plt.cm.tab20(np.linspace(0, 19, 20,dtype=int))
        colors2 = plt.cm.tab20b(np.linspace(0, 19, 20,dtype=int))

        # combine them and build a new colormap
        colors = np.vstack((colors1, colors2))
        cm = mcolors.ListedColormap(colors)

        

        for i1 in range(len(self.y_)):
            studyID_plt = int(study_ID[self.mask_meas_use][i1]) 
            if not_seen[studyID_plt]:
                axes[index_ax].plot(self.x_[i1],self.y_[i1],'o', color=cm(studyID_plt), label=df_keys[studyID_plt] )
                not_seen[studyID_plt] = False
            else:
                axes[index_ax].plot(self.x_[i1],self.y_[i1],'o', color=cm(studyID_plt) )
        
        if self.xaxis == 'modelled':
            axes[index_ax].set_ylabel(r'Measured [%s]' % units_string,fontsize=12)
            axes[index_ax].set_xlabel('Modelled [%s]' % units_string,fontsize=12) 
        else:
            axes[index_ax].set_xlabel(r'Measured [%s]' % units_string,fontsize=12)
            axes[index_ax].set_ylabel('Modelled [%s]' % units_string,fontsize=12) 
        
        try:
            minval = 0.9*min(min(self.x_[~(self.x_ == 0)]),min(self.y_[~(self.y_ == 0)]))
        except:
            minval = 1
            
        try:
            maxval = 1.11*max(max(self.x_),max(self.y_))
        except:
            maxval = 1e7
            
        axes[index_ax].plot([minval,maxval],[minval,maxval],'k-',label='1:1')
        
        
        if self.draw_errorbars:
            n_std = 2
            if units_ == '#/km2':
                log10_std_y = std_array[0]
            elif units_ == 'g/km2':
                log10_std_y = std_array[1]
                
            dy = n_std*log10_std_y
            
            y1u = minval+(minval*10**dy-minval)
            y2u = maxval+(maxval*10**dy-maxval)
            
            y1l = minval-(minval-10**(-dy)*minval)
            y2l = maxval-(maxval-10**(-dy)*maxval)
            
            axes[index_ax].plot([minval,maxval],[y1u,y2u],'r-',label='2x std from variogram')
            axes[index_ax].plot([minval,maxval],[y1l,y2l],'r-')
            
            # calculate % within errorbars
            mat = np.array([[minval,1],[maxval,1]])
            a_u,b_u = np.dot( np.linalg.inv(mat),np.array([[y1u],[y2u]]) )
            a_l,b_l = np.dot( np.linalg.inv(mat),np.array([[y1l],[y2l]]) )
            
            data_upper = a_u[0]*self.x_ + b_u[0]
            data_lower = a_l[0]*self.x_ + b_l[0]

            data_within = (self.y_ < data_upper) & (self.y_ > data_lower)
            percentage_within = 100*(np.sum(data_within)/len(data_within))

            axes[index_ax].set_title('R: %f \n percentage within %i$\sigma$: %3.1f' % (self.correlation,n_std,percentage_within) ,fontsize=12)
            
        axes[index_ax].set_yscale('log')
        axes[index_ax].set_xscale('log')
            
        
        return fig,axes 

#%%

def param_to_use(param_use,kde_use,units):
    """
    use_param = np.array([True,True,True,True]) #Use: beaching, sinking, source, river uncertainty

    t_beach_tau,t_sink_tau,t_sink_rate,sink_init,ratio_popToRiver,ratio_tourToRiver,ratio_fishToRiver,kde_kernel_bw,t_fragm_tau, = x_total
    """
    x_use = np.zeros(x0.shape,dtype=bool)
    if param_use[0]: #beaching
        x_use[0] = True
    if param_use[1]: #sinking
        x_use[1] = True
        x_use[2] = True
        x_use[3] = True
    if param_use[2]: #sources
        x_use[4] = True
        x_use[5] = True
    if param_use[3]: #use river uncertainty estimates
        x_use[7] = True
    if kde_use:
        x_use[6] = True

    return x_use


PLOT = False

projectFolder = '01_optimize/'
#dataDir = 'Data/Temp/34_noDelete_gemini/K100'
dataDir = 'Data/Temp/34_delete180Days_gemini/K10'

optimFile = 'optimizationParameters_180Days_K10'
progressFile = 'progressFile_180Days_K10'
final_results_file = 'Results_180Days_K10.pickle'

csvFiles = 'Data/Temp/34_CSVFiles_delete180Days_popRivFis_K10/'
#csvFiles = 'Data/Temp/34_CSVFiles_noDelete_popRivFis_K100/'

csvFiles_measurements = 'Data/Temp/34_CSVFiles_measurements/'

n_folders_per_type = 6
units = ['#/km2','g/km2']
std_array = [np.sqrt(0.1376), np.sqrt(0.2201)]


# parameters: beach_tau, sink_tau, sink_rate, sink_init, popRiver,fishRiver,KDE,river_factor
x_lower = np.array([ 1., 0.3,  3., 0.17, -1.3, -1.3,  0.05,  -1.]) #3* sigma
x_upper = np.array([ 3.,  1.72, 15., 0.44,  1.3,  1.3,  0.20, 1.])

log_scale = np.array([True,True,False,False,True,True,False,False]) 
                             
use_param = np.array([True,True,True,True]) #Use: beaching, sinking, source, degradation_tau, degradation_musig, river uncertainty, frag_dim
use_kde = True # optimize KDE as well
x_standard = np.array([10**2.5,10**0.3,2.00,0.0,1,1,0.2,-1.]) #standard parameters to use when on of 'use_param' is set to false. 


x0 = 0.5*np.ones(len(x_standard)) # init. optimization point 

CMAES_bounds = [0,1]

WRITE = True

use_SSD = True

time_origin = datetime(2006,2,1,0,0)


if os.environ['USER'] == 'mikael': # hp laptop
    import matplotlib.pyplot as plt

    homeDir = '/Volumes/externe_SSD/kaandorp/'

    outDir = os.path.join(os.getcwd(), projectFolder)
    
    particleFolders_pop = os.path.join(homeDir , dataDir, 'popOnly/')
    folders_pop = np.sort(glob.glob(os.path.join(particleFolders_pop,'*')))[0:n_folders_per_type]

    particleFolders_rivers = os.path.join(homeDir , dataDir, 'riversOnly/')
    folders_rivers = np.sort(glob.glob(os.path.join(particleFolders_rivers,'*')))[0:n_folders_per_type]    
    
    particleFolders_fisheries = os.path.join(homeDir , dataDir, 'fisheriesOnly/')
    folders_fisheries = np.sort(glob.glob(os.path.join(particleFolders_fisheries,'*')))[0:n_folders_per_type]   
    
    csvFiles = os.path.join(homeDir,csvFiles)
    csvFiles_measurements = os.path.join(homeDir,csvFiles_measurements)

if os.environ['USER'] == 'kaandorp': # desktop
    import matplotlib.pyplot as plt

    if use_SSD:
        homeDir = '/Volumes/externe_SSD/kaandorp/'
    else:
        homeDir = os.environ['HOME']
    outDir = os.path.join(os.getcwd(), projectFolder)
#    dataFile = os.path.join(homeDir, 'Documents', 'plastic_data_v3.xlsx')
    
    particleFolders_pop = os.path.join(homeDir , dataDir, 'popOnly/')
    folders_pop = np.sort(glob.glob(os.path.join(particleFolders_pop,'*')))[0:n_folders_per_type]

    particleFolders_rivers = os.path.join(homeDir , dataDir, 'riversOnly/')
    folders_rivers = np.sort(glob.glob(os.path.join(particleFolders_rivers,'*')))[0:n_folders_per_type]    
    
    particleFolders_fisheries = os.path.join(homeDir , dataDir, 'fisheriesOnly/')
    folders_fisheries = np.sort(glob.glob(os.path.join(particleFolders_fisheries,'*')))[0:n_folders_per_type]   
    
    csvFiles = os.path.join(homeDir,csvFiles)
    csvFiles_measurements = os.path.join(homeDir,csvFiles_measurements)

if os.environ['USER'] == 'kaand004': #gemini cluster

    os.nice(0) #set low priority for cpu use
    os.environ['MKL_NUM_THREADS'] = '4'
    
    homeDir = '/scratch/kaand004'
    outDir = os.path.join(homeDir, projectFolder)
    
    particleFolders_pop = os.path.join(homeDir , dataDir, 'popOnly/')
    folders_pop = np.sort(glob.glob(os.path.join(particleFolders_pop,'*')))[0:n_folders_per_type]
    
    particleFolders_rivers = os.path.join(homeDir , dataDir, 'riversOnly/')
    folders_rivers = np.sort(glob.glob(os.path.join(particleFolders_rivers,'*')))[0:n_folders_per_type]    
    
    particleFolders_fisheries = os.path.join(homeDir , dataDir, 'fisheriesOnly/')
    folders_fisheries = np.sort(glob.glob(os.path.join(particleFolders_fisheries,'*')))[0:n_folders_per_type]   
    
    csvFiles = os.path.join(homeDir,csvFiles)
    csvFiles_measurements = os.path.join(homeDir,csvFiles_measurements)


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
init_source[index_plus[len(folders_pop)+len(folders_rivers)]:] = 3

lons = np.linspace(-6,36.25,677)
lats = np.linspace(30.1875,45.9375,253)
init_lon_id = np.digitize(init_lon,lons,right=True)
init_lat_id = np.digitize(init_lat,lats,right=True)

time_20100101 = (datetime(2010,1,1,0,0)-time_origin).total_seconds()
time_20101231 = (datetime(2010,12,31,23,59)-time_origin).total_seconds()

n_pop = len(init_source[init_source == 1])
n_river = len(init_source[init_source == 2])
n_fish = len(init_source[init_source == 3])

n_river_2010 = len(init_source[(init_source == 2) & (init_time > time_20100101) & (init_time < time_20101231)])
kg_year_rivers_mid = 1248000 # from Lebreton (2017), mid estimate
kg_particle_river_mid = kg_year_rivers_mid / n_river_2010 

init_month_id = np.array([ (time_origin + timedelta(seconds=int(init_time[i1]))).month - 1 for i1 in range(len(init_time.values))])
riverInput_low = np.loadtxt('datafile_riverInputMatrix_677x_253y_low').reshape([12,len(lats),len(lons)])
riverInput_mid = np.loadtxt('datafile_riverInputMatrix_677x_253y').reshape([12,len(lats),len(lons)])
riverInput_high = np.loadtxt('datafile_riverInputMatrix_677x_253y_high').reshape([12,len(lats),len(lons)])


filenames = np.sort(glob.glob(os.path.join(csvFiles,'*particles.csv')))

# cast to appropriate types, check:
"""
int_types = ["uint8", "int8", "int32"]
for it in int_types:
    print(np.iinfo(it))
"""
dtype_pd = dtype={'# ID': np.int32,'lon':np.float32,'lat':np.float32,'coastalAge_1cell':np.float32,'coastalAge_2cell':np.float32,
                  'age':np.float32,'dot_curr':np.float32,'dot_wind':np.float32,'dot_wind_Ek':np.float32,'dot_stokes':np.float32,
                  'mag_curr':np.float32,'mag_stokes':np.float32,'mag_wind':np.float32}
usecols_pd = usecols=['# ID', 'lon', 'lat', 'coastalAge_1cell', 'age']

dfs = [delayed(load)(fn) for fn in filenames]

df = dd.from_delayed(dfs).compute()


# in the entire dataframe, save starting indices corresponding to each measurement
indices_0 = np.where(df.index == 0)[0]

indices_0 = np.append(indices_0,df.shape[0])

x_use = param_to_use(use_param,use_kde,units)

popDen_max = init_popDen.max()

measured, lon_meas, lat_meas, study_ID, basin_ID, kukulka_meas, datetime_meas, unit_meas, sizeRangeCorr_meas = initMeasured()


#%% 

def calc_dgdm(m,eps,type_='fwd'):
    """
    function calculating the Jacobian of the forward model using finite differences
    
    m: array with parameters
    eps: finite difference step
    type_: finite difference type (fwd or central)
    """
    
    
    def dgdmi_finite_diff(i_m,m,eps,type_='central',g_m=None):
        mp = m.copy()
        mp[i_m] += eps
        mm = m.copy()
        mm[i_m] -= eps
        
        if type_ == 'central': 
            print(mp)
            print(mm)
            
            g_mp,_,_,_,_,success_p,_,_,_ = forward_model(mp).calculate_matrices()
            g_mm,_,_,_,_,success_m,_,_,_ = forward_model(mm).calculate_matrices()
        
            dg_dm_i = (g_mp - g_mm) / (2*eps)
        elif type_ == 'fwd':
            g_mp,_,_,_,_,success_p,_,_,_ = forward_model(mp).calculate_matrices()
            dg_dm_i = (g_mp - g_m) / eps
        else:
            raise RuntimeError('unknown fd type')
    
        return dg_dm_i
    
    if type_ == 'fwd':
        g_m,_,_,_,_,success_g,_,_,_ = forward_model(m).calculate_matrices()
    else:
        g_m = None
    
    for i1 in range(len(x0)):
            
        if i1 == 0:
            dgdm_mat = dgdmi_finite_diff(0,m,eps,type_,g_m)
        else:
            if x_use[i1] == True:
                dgdm_mat = np.vstack((dgdm_mat,dgdmi_finite_diff(i1,m,eps,type_,g_m)))
            else: #not use this param: set gradient to zero
                if i1 == 1:
                    n_ = len(dgdm_mat)
                else:
                    n_ = dgdm_mat.shape[1]
                dgdm_mat = np.vstack((dgdm_mat,np.zeros(n_)))

        print('%i / %i' % (i1,len(x0)))
        
        f = open(os.path.join(outDir,progressFile), 'a')
        f.write('%i / %i\n' % (i1,len(x0)))
        f.close() 
        
    return dgdm_mat.T
    
#%% Quasi-newton optimization

param_names = ['tau_beach','tau_sink','rate_sink','sink_init','pop:riv','fis:riv','KDE','tau_frag.','mu_frag.','sig_frag.','f_river','d_frag.']

m_prior = 0.5*np.ones(x0.shape)


x_standard_unitless = forward_model(x0).to_unitless(x_standard,x_lower,x_upper,log_scale)
m_prior[~x_use] = x_standard_unitless[~x_use] # fix fragmentation_dimension

#m_0 = 0.5*np.ones(x0.shape)
#m_0[~x_use] = x_standard_unitless[~x_use]

#previously found starting point to speed up computations
m_0 = np.array([ 0.26098583,  0.24523944, 0.3,  0.5,  0.66775254,        0.36910061,  1.06042558,  0.     ])

n_iter = 20
eps = 0.005
mu = .8
mu_end = 0.1
m_n = m_0.copy()

f_mu = (mu_end/mu)**(1/n_iter)

for i1 in range(n_iter):
    dgdm_mat = calc_dgdm(m_n,eps,type_='fwd')
    g_m,d,C_d,C_m,param_deviation,all_success,_,_,_ = forward_model(m_n).calculate_matrices()
    t1 = np.linalg.inv(np.dot(dgdm_mat.T,np.dot(np.linalg.inv(C_d),dgdm_mat)) + np.linalg.inv(C_m))
    t2 = np.dot(dgdm_mat.T,np.dot(np.linalg.inv(C_d),(g_m-d))) 
    t3 = np.dot(np.linalg.inv(C_m),(m_n - m_prior))
    m_np1 = m_n - mu*np.dot(t1,(t2+t3))
    
    m_unitfull= forward_model(x0).to_unitfull(m_np1,x_lower,x_upper,log_scale)
    print('%7.7s|%7.7s|%7.7s|%7.7s|%7.7s|%7.7s|%7.7s|%7.7s|%7.7s|%7.7s|%7.7s|%7.7s' % tuple(param_names)   ) 
    print('%4.4s[d]|%4.4s[w]|%4.4s[-]|%4.4s[-]|%4.4s   |%4.4s   |%4.4s   |%4.4s[d]|%4.4s   |%4.4s   |%4.4s   |%4.4s   ' % tuple(m_unitfull)   )     
    print(np.array_repr(m_np1).replace('\n', ''))
    
    f = open(os.path.join(outDir,progressFile), 'a')
    f.write('%7.7s|%7.7s|%7.7s|%7.7s|%7.7s|%7.7s|%7.7s|%7.7s|%7.7s|%7.7s|%7.7s|%7.7s\n' % tuple(param_names)   ) 
    f.write('%4.4s[d]|%4.4s[w]|%4.4s[-]|%4.4s[-]|%4.4s   |%4.4s   |%4.4s   |%4.4s[d]|%4.4s   |%4.4s   |%4.4s   |%4.4s   \n' % tuple(m_unitfull)   )    
    f.write(np.array_repr(m_np1).replace('\n', '')+'\n')
    f.close()  
    
    m_n = m_np1
     
    #update relaxation factor
    mu = mu*f_mu
    
        


#%% Plot results
    
vals_labels = [r'$\tau_{beach}$ [days]',r'$\tau_{sink}$ [days]',r'$r_{sink}$ [-]',r'sink$_{init.}$ [-]','pop:river [-]','fish:river [-]','KDE$_{bw}$ [-]',r'$\tau_{fragm.}$ [days]',r'$\mu_{fragm.}$ [ln(mm)]',r'$\sigma_{fragm.}$ [ln(mm)]','river$_{low-high}$ [-]','$d_{fragm.}$ [-]']
to_plot = [True,True,True,True,True,True,True,False,False,False,True,False]

eps = 0.005
m_prior = 0.5*np.ones(x0.shape)
x_op_dimless = m_n.copy()
g_m,d,C_d,C_m,param_deviation,all_success,mask_measurements,measurements,modelled = forward_model(x_op_dimless).calculate_matrices()    
dgdm_mat = calc_dgdm(x_op_dimless,eps,type_='fwd')

C_m_posterior = np.linalg.inv(np.dot(dgdm_mat.T,np.dot(np.linalg.inv(C_d),dgdm_mat)) + np.linalg.inv(C_m))

sigma_posterior = np.sqrt(np.diag(C_m_posterior))

mismatch_2S = np.dot((g_m-d).T,np.dot(np.linalg.inv(C_d),(g_m-d))) + np.dot((m_n-m_prior)[x_use].T,np.dot(np.linalg.inv((np.eye(x_use.sum())*C_m.diagonal()[x_use])),(m_n-m_prior)[x_use]))

dimless_low = np.floor(x_op_dimless - 3*sigma_posterior)
dimless_high = np.ceil(x_op_dimless + 3*sigma_posterior)

if PLOT:
    n_rows = int(np.ceil(sum(to_plot)/2))
    fig,ax = plt.subplots(n_rows,2,figsize=(15,15))
    for i1 in range(sum(to_plot)):
        col_ = i1 % 2
        row_ = i1 // 2
        
    
        x_ = np.linspace(dimless_low[to_plot][i1],dimless_high[to_plot][i1],1000)
            
        x_scaled = x_lower[to_plot][i1]+x_*(x_upper[to_plot][i1]-x_lower[to_plot][i1])
        if log_scale[to_plot][i1]:
            x_scaled = 10**x_scaled
        
        if i1 == 1: #weeks to days
            x_scaled *= 7
        
        
        pdf_ = norm(x_op_dimless[to_plot][i1],sigma_posterior[to_plot][i1])
        y_ = pdf_.pdf(x_)
        
        if log_scale[to_plot][i1]:
            ax[row_,col_].semilogx(x_scaled,y_,'k-')
        else:
            ax[row_,col_].plot(x_scaled,y_,'k-')
        
        
        ax[row_,col_].text(0.7,0.8,np.array(vals_labels)[to_plot][i1],transform=ax[row_,col_].transAxes,fontsize=14)

#%% Save results
        
dict_results = {}
dict_results['x_op_dimless'] = m_n.copy()
dict_results['sigma_posterior'] = sigma_posterior.copy()
dict_results['dgdm_mat'] = dgdm_mat.copy()
dict_results['C_d'] = C_d.copy()
dict_results['C_m'] = C_m.copy()
dict_results['x_upper'] = x_upper.copy()
dict_results['x_lower'] = x_lower.copy()
dict_results['log_scale'] = log_scale.copy()
dict_results['mismatch_2S'] = mismatch_2S.copy()
dict_results['mask_measurements'] = mask_measurements.copy()
dict_results['measurements'] = measurements.copy()
dict_results['modelled'] = modelled.copy()
outfile = open(os.path.join(outDir,final_results_file),'wb')
pickle.dump(dict_results,outfile)
outfile.close()

dict_results_2 = {}
dict_results_2['C_m_post'] = C_m_posterior.copy()
dict_results_2['x_op_dimless'] = m_n.copy()