# -*- coding: utf-8 -*-
"""
======
HBV-96d
======

ditributed-Lumped hydrological model.

This is the HBV-96 implementation by Juan Chacon at IHE-Delft, NL. This code
implements a distributed-lumped version of the HBV-96, as described in 
Chacon-Hurtado, et al. (2017).

The model uses the path length to route the superficial water, while the lower 
zone is assumed as lumped.

@author: Juan Carlos Chacon-Hurtado (jc.chaconh@gmail.com)                                  

Version
-------
13-11-2017 - V_0.0 - First implementation
"""
from __future__ import division, print_function
import numpy as np
import scipy.optimize as opt
import gdal

# HBV base model parameters
P_LB = [-1.5, #ltt
        0.001, #utt
        0.001, #ttm
        0.04, #cfmax [mm c^-1 h^-1]
        50.0, #fc
        0.6, #ecorr
        0.001, #etf
        0.2, #lp
        0.00042, #k [h^-1] upper zone
        0.0000042, #k1 lower zone
        0.001, #alpha
        1.0, #beta
        0.001, #cwh
        0.01, #cfr
        0.0, #c_flux
        0.001, #perc mm/h
        0.6, #rfcf
        0.4,  #sfcf
        1] # Maxbas

P_UB = [2.5, #ttm
        3.0, #utt
        2.0, #ttm
        0.4, #cfmax [mm c^-1 h^-1]
        500.0, #fc
        1.4, #ecorr
        5.0, #etf
        0.5, #lp
        0.0167, #k upper zone
        0.00062, #k1 lower zone
        1.0, #alpha
        6.0, #beta
        0.1, #cwh
        1.0, #cfr
        0.08, #c_flux - 2mm/day
        0.125, #perc mm/hr
        1.4, #rfcf
        1.4, #sfcf
        10] # maxbas

DEF_ST = [0.0, 10.0, 10.0, 10.0, 0.0]
DEF_q0 = 10.0

# Get random parameter set
def get_random_pars():
    return np.random.uniform(P_LB, P_UB)

def _precipitation(temp, ltt, utt, prec, rfcf, sfcf, tfac):
    '''
    ==============
    Precipitation
    ==============

    Precipitaiton routine of the HBV96 model.

    If temperature is lower than ltt, all the precipitation is considered as
    snow. If the temperature is higher than utt, all the precipitation is
    considered as rainfall. In case that the temperature is between ltt and
    utt, precipitation is a linear mix of rainfall and snowfall.

    Parameters
    ----------
    temp : float
        Measured temperature [C]
    ltt : float
        Lower temperature treshold [C]
    utt : float
        Upper temperature treshold [C]
    prec : float 
        Precipitation [mm]
    rfcf : float
        Rainfall corrector factor
    sfcf : float
        Snowfall corrector factor

    Returns
    -------
    _rf : float
        Rainfall [mm]
    _sf : float
        Snowfall [mm]
    '''

    if temp <= ltt:
        _rf = 0.0
        _sf = prec*sfcf

    elif temp >= utt:
        _rf = prec*rfcf
        _sf = 0.0

    else:
        _rf = ((temp-ltt)/(utt-ltt)) * prec * rfcf
        _sf = (1.0-((temp-ltt)/(utt-ltt))) * prec * sfcf

    return _rf, _sf


def _snow(cfmax, tfac, temp, ttm, cfr, cwh, _rf, _sf, wc_old, sp_old):
    '''
    ====
    Snow
    ====
    
    Snow routine of the HBV-96 model.
    
    The snow pack consists of two states: Water Content (wc) and Snow Pack 
    (sp). The first corresponds to the liquid part of the water in the snow,
    while the latter corresponds to the solid part. If the temperature is 
    higher than the melting point, the snow pack will melt and the solid snow
    will become liquid. In the opposite case, the liquid part of the snow will
    refreeze, and turn into solid. The water that cannot be stored by the solid
    part of the snow pack will drain into the soil as part of infiltration.

    Parameters
    ----------
    cfmax : float 
        Day degree factor
    tfac : float
        Temperature correction factor
    temp : float 
        Temperature [C]
    ttm : float 
        Temperature treshold for Melting [C]
    cfr : float 
        Refreezing factor
    cwh : float 
        Capacity for water holding in snow pack
    _rf : float 
        Rainfall [mm]
    _sf : float 
        Snowfall [mm]
    wc_old : float 
        Water content in previous state [mm]
    sp_old : float 
        snow pack in previous state [mm]

    Returns
    -------
    _in : float 
        Infiltration [mm]
    _wc_new : float 
        Water content in posterior state [mm]
    _sp_new : float 
        Snowpack in posterior state [mm]
    '''

    if temp > ttm:

        if cfmax*(temp-ttm) < sp_old+_sf:
            _melt = cfmax*(temp-ttm)
        else:
            _melt = sp_old+_sf

        _sp_new = sp_old + _sf - _melt
        _wc_int = wc_old + _melt + _rf

    else:
        if cfr*cfmax*(ttm-temp) < wc_old+_rf:
            _refr = cfr*cfmax*(ttm - temp)
        else:
            _refr = wc_old + _rf

        _sp_new = sp_old + _sf + _refr
        _wc_int = wc_old - _refr + _rf

    if _wc_int > cwh*_sp_new:
        _in = _wc_int-cwh*_sp_new
        _wc_new = cwh*_sp_new
    else:
        _in = 0.0
        _wc_new = _wc_int

    return _in, _wc_new, _sp_new


def _soil(fc, beta, etf, temp, tm, e_corr, lp, tfac, c_flux, inf,
          ep, sm_old, uz_old):
    '''
    ====
    Soil
    ====
    
    Soil routine of the HBV-96 model.
    
    The model checks for the amount of water that can infiltrate the soil, 
    coming from the liquid precipitation and the snow pack melting. A part of 
    the water will be stored as soil moisture, while other will become runoff, 
    and routed to the upper zone tank.

    Parameters
    ----------
    fc : float 
        Filed capacity
    beta : float 
        Shape coefficient for effective precipitation separation
    etf : float 
        Total potential evapotranspiration
    temp : float 
        Temperature
    tm : float 
        Average long term temperature
    e_corr : float 
        Evapotranspiration corrector factor
    lp : float _soil 
        wilting point
    tfac : float 
        Time conversion factor
    c_flux : float 
        Capilar flux in the root zone
    _in : float 
        actual infiltration
    ep : float 
        actual evapotranspiration
    sm_old : float 
        Previous soil moisture value
    uz_old : float 
        Previous Upper zone value

    Returns
    -------
    sm_new : float 
        New value of soil moisture
    uz_int_1 : float 
        New value of direct runoff into upper zone
    '''

    qdr = max(sm_old + inf - fc, 0)
    _in = inf - qdr
    _r = ((sm_old/fc)** beta) * _in
    _ep_int = max((1.0 + etf*((temp - tm)/tm))*e_corr*ep, 0)
    _ea = min(_ep_int, (sm_old/(lp*fc))*_ep_int)

    _cf = c_flux*((fc - sm_old)/fc)
    sm_new = max(sm_old + _in - _r + _cf - _ea, 0)
    uz_int_1 = uz_old + _r - _cf + qdr

    return sm_new, uz_int_1


def _response(tfac, perc, alpha, k, k1, area, lz_old, uz_int_1):
    '''
    ========
    Response
    ========
    The response routine of the HBV-96 model.
    
    The response routine is in charge of transforming the current values of 
    upper and lower zone into discharge. This routine also controls the 
    recharge of the lower zone tank (baseflow). The transformation of units 
    also occurs in this point.
    
    Parameters
    ----------
    tfac : float
        Number of hours in the time step
    perc : float
        Percolation value [mm\hr]
    alpha : float
        Response box parameter
    k : float
        Upper zone recession coefficient
    k1 : float 
        Lower zone recession coefficient
    area : float
        Catchment area [Km2]
    lz_old : float 
        Previous lower zone value [mm]
    uz_int_1 : float 
        Previous upper zone value before percolation [mm]
    qdr : float
        Direct runoff [mm]
    
    '''    
    lz_int_1 = lz_old + np.min([perc, uz_int_1])
    uz_int_2 = np.max([uz_int_1 - perc, 0.0])

    _q_0 = k*(uz_int_2**(1.0 + alpha))
    _q_1 = k1*lz_int_1

    uz_new = max(uz_int_2 - (_q_0), 0)
    lz_new = max(lz_int_1 - (_q_1), 0)

    q_new = area*(_q_0 + _q_1)/(3.6*tfac)

    return q_new, uz_new, lz_new, uz_int_2, lz_int_1


def _tf(maxbas):
    ''' Transfer function weight generator '''
    wi = []
    for x in range(1, maxbas+1):
        if x <= (maxbas)/2.0:
            # Growing transfer
            wi.append((x)/(maxbas+2.0))
        else:
            # Receding transfer
            wi.append(1.0 - (x+1)/(maxbas+2.0))
    
    #Normalise weights
    wi = np.array(wi)/np.sum(wi)
    return wi

def _routing(q, maxbas=1):
    """
    This function implements the transfer function using a triangular 
    function
    """
    assert maxbas >= 1, 'Maxbas value has to be larger than 1'
    # Get integer part of maxbas
    maxbas = int(maxbas)
    
    # get the weights
    w = _tf(maxbas)
    
    # rout the discharge signal
    q_r = np.zeros_like(q, dtype='float64')
    q_temp = q
    for w_i in w:
        q_r += q_temp*w_i
        q_temp = np.insert(q_temp, 0, 0.0)[:-1]

    return q_r


def _step_run(p, p2, v, St, extra_out=False):
    '''
    ========
    Step run
    ========
    
    Makes the calculation of next step of discharge and states
    
    Parameters
    ----------
    p : array_like [18]
        Parameter vector, set up as:
        [ltt, utt, ttm, cfmax, fc, ecorr, etf, lp, k, k1, 
        alpha, beta, cwh, cfr, c_flux, perc, rfcf, sfcf]
    p2 : array_like [2]
        Problem parameter vector setup as:
        [tfac, area]
    v : array_like [4]
        Input vector setup as:
        [prec, temp, evap, llt]
    St : array_like [5]
        Previous model states setup as:
        [sp, sm, uz, lz, wc]

    Returns
    -------
    q_new : float
        Discharge [m3/s]
    St : array_like [5]
        Posterior model states
    '''
    ## Parse of parameters from input vector to model
    ltt = p[0]
    utt = p[1]
    ttm = p[2]
    cfmax = p[3]
    fc = p[4]
    e_corr = p[5]
    etf = p[6]
    lp = p[7]
    k = p[8]
    k1 = p[9]
    alpha = p[10]
    beta = p[11]
    cwh = p[12]
    cfr = p[13]
    c_flux = p[14]
    perc = p[15]
    rfcf = p[16]
    sfcf = p[17]
    
    ## Non optimisable parameters
    tfac = p2[0]
    area = p2[1]

    ## Parse of Inputs
    avg_prec = v[0] # Precipitation [mm]
    temp = v[1] # Temperature [C]
    ep = v[2] # Long terms (monthly) Evapotranspiration [mm]
    tm = v[3] #Long term (monthly) average temperature [C]

    ## Parse of states
    sp_old = St[0]
    sm_old = St[1]
    uz_old = St[2]
    lz_old = St[3]
    wc_old = St[4]

    rf, sf = _precipitation(temp, ltt, utt, avg_prec, rfcf, sfcf, tfac)
    inf, wc_new, sp_new = _snow(cfmax, tfac, temp, ttm, cfr, cwh, rf, sf,
                               wc_old, sp_old)
    sm_new, uz_int_1 = _soil(fc, beta, etf, temp, tm, e_corr, lp, tfac, c_flux,
                            inf, ep, sm_old, uz_old)
    q_new, uz_new, lz_new, uz_int_2, lz_int_1 = _response(tfac, perc, alpha, 
                                                          k, k1, area, lz_old,
                                                          uz_int_1)
    
    return q_new, [sp_new, sm_new, uz_new, lz_new, wc_new], uz_int_2, lz_int_1


def simulate(avg_prec, temp, et, par, p2, init_st=None, ll_temp=None, 
             q_0=DEF_q0, extra_out=False):
    '''
    ========
    Simulate
    ========

    Run the HBV model for the number of steps (n) in precipitation. The
    resluts are (n+1) simulation of discharge as the model calculates step n+1

    
    Parameters
    ----------
    avg_prec : array_like [n]
        Average precipitation [mm/h]
    temp : array_like [n]
        Average temperature [C]
    et : array_like [n]
        Potential Evapotranspiration [mm/h]
    par : array_like [18]
        Parameter vector, set up as:
        [ltt, utt, ttm, cfmax, fc, ecorr, etf, lp, k, k1, 
        alpha, beta, cwh, cfr, c_flux, perc, rfcf, sfcf]
    p2 : array_like [2]
        Problem parameter vector setup as:
        [tfac, area]
    init_st : array_like [5], optional
        Initial model states, [sp, sm, uz, lz, wc]. If unspecified, 
        [0.0, 30.0, 30.0, 30.0, 0.0] mm
    ll_temp : array_like [n], optional
        Long term average temptearature. If unspecified, calculated from temp.
    q_0 : float, optional
        Initial discharge value. If unspecified set to 10.0
    

    Returns
    -------
    q_sim : array_like [n]
        Discharge for the n time steps of the precipitation vector [m3/s]
    st : array_like [n, 5]
        Model states for the complete time series [mm]
    '''

    if init_st is None:
        st = [DEF_ST, ]

    if ll_temp is None:
        ll_temp = [np.mean(temp), ] * len(avg_prec)

    q_sim = [q_0, ]
    
    #print(st)
    uz_int_2 = [st[0][2], ]
    lz_int_1 = [st[0][3], ]
    
    for i in xrange(len(avg_prec)):
        v = [avg_prec[i], temp[i], et[i], ll_temp[i]]
        q_out, st_out, uz_int_2_out, lz_int_1_out = _step_run(par, p2, v, st[i])
        q_sim.append(q_out)
        st.append(st_out)
        uz_int_2.append(uz_int_2_out)
        lz_int_1.append(lz_int_1_out)
    
    if len(p2) > 2:  # Forcing maxbas to be predefined
        maxbas = p2[2]  
    elif len(par) > 18:  # Putting maxbas as parameter to be optimised
        maxbas = par[18]
    else:
        maxbas = 1
    
    q_tr = _routing(np.array(q_sim), maxbas)
    
    if extra_out:
        return q_tr, st, uz_int_2, lz_int_1
    
    else:
        return q_tr, st

def _add_mask(var, dem=None, mask=None, no_val=None):
    '''
    Put a mask in the spatially distributed values
    
    Parameters
    ----------
    var : nd_array
        Matrix with values to be masked
    cut_dem : gdal_dataset
        Instance of the gdal raster of the catchment to be cutted with. DEM 
        overrides the mask_vals and no_val
    mask_vals : nd_array
        Mask with the no_val data
    no_val : float
        value to be defined as no_val. Will mask anything is not this value
    
    Returns
    -------
    var : nd_array
        Array with masked values 
    '''
    
    if dem is not None:
        mask, no_val = _get_mask(dem)
    
    # Replace the no_data value
    assert var.shape == mask.shape, 'Mask and data do not have the same shape'
    var[mask == no_val] = no_val
    
    return var

def _get_mask(dem):
    no_val = dem.GetRasterBand(1).GetNoDataValue()
    mask = dem.ReadAsArray()
    return mask, no_val

def _get_targets(dem):
    '''
    Returns the centres of the interpolation targets
    
    Parameters
    ----------
    dem : gdal_Dataset
        Get the data from the gdal datasetof the DEM
    
    Returns
    -------
    
    coords : nd_array [nxm - nan, 2]
        Array with a list of the coordinates to be interpolated, without the Nan
    
    mat_range : nd_array [n, m]
        Array with all the centres of cells in the domain of the DEM (rectangular)
    
    '''
    # Getting data for the whole grid
    x_init, xx_span, xy_span, y_init, yy_span, yx_span = dem.GetGeoTransform()
    shape_dem = dem.ReadAsArray().shape
    
    # Getting data of the mask
    no_val = dem.GetRasterBand(1).GetNoDataValue()
    mask = dem.ReadAsArray()
    
    # Adding 0.5 to get the centre
    x = np.array([x_init + xx_span*(i+0.5) for i in xrange(shape_dem[0])])
    y = np.array([y_init + yy_span*(i+0.5) for i in xrange(shape_dem[1])])
    #mat_range = np.array([[(xi, yi) for xi in x] for yi in y])
    mat_range = [[(xi, yi) for xi in x] for yi in y]
    
    # applying the mask
    coords = []
    for i in xrange(len(x)):
        for j in range(len(y)):
            if mask[j, i] != no_val:
                coords.append(mat_range[j][i])
                #mat_range[j, i, :] = [np.nan, np.nan]

    return np.array(coords), np.array(mat_range)


def distributed(flp, sp_prec, sp_et, sp_temp, sp_pars, p2, init_st=None, 
                ll_temp=None, q_0=None):
    '''
    Runs the Distributed-Lumped HBV96.
    
    Parameters
    ----------
    
    flp : GDAL instance 
        Raster shapefile with maximum flow path length
    sp_prec : nd_array [t, x, y]
        Array with precipitation data [t, x, y]. Has to be compatible with the 
        size of flp
    sp_et : nd_array [t, x, y]
        Array with evapotranspiration data [t, x, y]. Has to be compatible 
        with the size of flp
    sp_temp : nd_array [t, x, y]
        Array with temperature data [t, x, y]. Has to be compatible with the 
        size of flp
    sp_pars : nd_array [x, y, 17] or [17]
        Array with the spatial distribution of the Parameter set. If a single 
        vector is passed as parameter set, it will be assumed as constant
    p2 : nd_array [3]
        Array with the non-optimisable parameters, TFAC, AREA and (maximum) 
        MAXBAS
    init_st : nd_array [x, y, 5] - [optional]
        Array with the initial states of the system. If the array is 1d, the 
        initial states of the system are considered homogeneous.
    ll_temp : nd_array [t, x, y] - [optional]
        Array with the spatial distribution of the long-term average 
        temperature. Has to be compatible with the size of flp
    q_0 : nd_array [x,y] - [optional]
        Array with the spatial estimates of initial discharge values
        
    Return
    ------
    q_out : nd_array [t + 1]
        Vector with discharge at the outlet
    st : nd_array [t, x, y]
        Matrix with states of the 
    
    '''
    
    #%%
#    init_st = None
#    q_0 = None
#    sp_pars = P_UB
#    p2 = [1, 145.0]
    # initialise Q
    
    if q_0 is None:
        q_0 = DEF_q0
        
    n_steps = sp_prec.shape[0] + 1
    # intiialise vector of nans to fill states
    _dummy_states = np.empty([n_steps, 5])
    _dummy_states[:] = np.nan
    
    # Get the mask
    mask, no_val = _get_mask(flp)
    _x_ext, _y_ext = mask.shape
    
    # Parse the parameter to its spatially distributed form
    sp_pars = np.array(sp_pars)
    if len(sp_pars.shape) == 1:  # Checking for distributed pars
        
        # assign the maxbas value
        if len(p2) > 2:  # Forcing maxbas to be predefined
            maxbas = p2[2]  
        elif len(sp_pars) > 18:  # Putting maxbas as parameter to be optimised
            maxbas = sp_pars[18]
        else:
            maxbas = 1
        
        sp_pars = np.tile(sp_pars, [_x_ext, _y_ext, 1])
    
    else:
        if len(p2) > 2:  # Forcing maxbas to be predefined
            maxbas = p2[2]  
        elif len(sp_pars[0,0,:]) > 18:  # Putting maxbas as parameter to be optimised
            maxbas = sp_pars[18]
        else:
            maxbas = 1
    
    #maxbas = 4
    
    # apply routing to UZ. farthest location gets double of the maxbas, 
    # closest gets min (1)
    # get the distance vector from the flp UZ = 2
    dist_map = flp.ReadAsArray()
    dist_map[mask == no_val] = np.nan
    
    max_dist = np.nanmax(dist_map)
    min_dist = np.nanmin(dist_map)
    resize_fun = lambda (x): np.round(((((x - min_dist)/(max_dist - min_dist))*(1*maxbas - 1)) + 1), 0)
    
    # maxbas_map
    maxbas_map = resize_fun(dist_map)
    
    st = []  # Spatially distributed states
    q_lz = []
    q_uz = []
    
    # TODO Make function here to add to multiprocess
    for x in xrange(_x_ext):
        st_i = []
        q_lzi = []
        q_uzi = []
#        q_out_i = []
        for y in xrange(_y_ext):
            if mask [x, y] != no_val:  # only for cells in the domain
                # Calculate the states per cell
                # TODO optimise for multiprocessing these loops   
                
                _, _st, _uzg, _lzg = simulate(avg_prec = sp_prec[:, x, y], 
                                              temp = sp_temp[:, x, y], 
                                              et = sp_et[:, x, y], 
                                              par = sp_pars[x, y, :], 
                                              p2 = p2, 
                                              init_st = init_st, 
                                              ll_temp = None, 
                                              q_0 = q_0,
                                              extra_out = True)
#                
                st_i.append(np.array(_st))
                q_lzi.append(np.array(sp_pars[x, y, 9])*_lzg)
                q_uz_temp = np.array(sp_pars[x, y, 8])*(np.power(_uzg, (1.0 + sp_pars[x, y, 10])))
                q_uzi.append(_routing(q_uz_temp, maxbas_map[x,y]))

            else:
                # Fill the empty spaces with a nan vector
                st_i.append(_dummy_states)
                q_lzi.append(_dummy_states[:,0])
                q_uzi.append(_dummy_states[:,0])

        st.append(st_i)
        q_lz.append(np.array(q_lzi))
        q_uz.append(np.array(q_uzi))
            
    st = np.array(st)
    q_lz = np.array(q_lz)
    q_lz = np.array([np.nanmean(q_lz[:,:,i]) for i in xrange(n_steps)])
    
    q_uz = np.array(q_uz)
    q_uz = np.array([np.nanmean(q_uz[:,:,i]) for i in xrange(n_steps)])
    
    q_out = (q_lz + q_uz) * p2[1]/ (p2[0]*3.6)

    return q_out, st

def _nse(q_rec, q_sim):
    '''
    ===
    NSE
    ===
    
    Nash-Sutcliffe efficiency. Metric for the estimation of performance of the 
    hydrological model
    
    Parameters
    ----------
    q_rec : array_like [n]
        Measured discharge [m3/s]
    q_sim : array_like [n] 
        Simulated discharge [m3/s]
        
    Returns
    -------
    f : float
        NSE value
    '''
    a = np.square(np.subtract(q_rec, q_sim))
    b = np.square(np.subtract(q_rec, np.nanmean(q_rec)))
    if a.any < 0.0:
        return(np.nan)
    f = 1.0 - (np.nansum(a)/np.nansum(b))
    return f


def _rmse(q_rec,q_sim):
    '''
    ====
    RMSE
    ====
    
    Root Mean Squared Error. Metric for the estimation of performance of the 
    hydrological model.
    
    Parameters
    ----------
    q_rec : array_like [n]
        Measured discharge [m3/s]
    q_sim : array_like [n] 
        Simulated discharge [m3/s]
        
    Returns
    -------
    f : float
        RMSE value
    '''
    erro = np.square(np.subtract(q_rec,q_sim))
    if erro.any < 0:
        return(np.nan)
    f = np.sqrt(np.nanmean(erro))
    return f

def calibrate(flow, avg_prec, temp, et, p2, init_st=None, ll_temp=None,
              x_0=None, x_lb=P_LB, x_ub=P_UB, obj_fun=_rmse, wu=10,
              verbose=False, tol=0.001, minimise=True, fun_nam='RMSE'):
    '''
    ## Not tested/implemented yet!
    =========
    Calibrate
    =========

    Run the calibration of the HBV-96. The calibration is used to estimate the
    optimal set of parameters that minimises the difference between 
    observations and modelled discharge.
    
    Parameters
    ----------
    
    flow : array_like [n]
        Measurements of discharge [m3/s]
    avg_prec : array_like [n]
        Average precipitation [mm/h]
    temp : array_like [n]
        Average temperature [C]
    et : array_like [n]
        Potential Evapotranspiration [mm/h] 
    p2 : array_like [2]
        Problem parameter vector setup as:
        [tfac, area]
    init_st : array_like [5], optional
        Initial model states, [sp, sm, uz, lz, wc]. If unspecified, 
        [0.0, 30.0, 30.0, 30.0, 0.0] mm
    ll_temp : array_like [n], optional
        Long term average temptearature. If unspecified, calculated from temp.
    x_0 : array_like [18], optional
        First guess of the parameter vector. If unspecified, a random value
        will be sampled between the boundaries of the 
    x_lb : array_like [18], optional
        Lower boundary of the parameter vector. If unspecified, a random value
        will be sampled between the boundaries of the 
    x_ub : array_like [18], optional
        First guess of the parameter vector. If unspecified, a random value
        will be sampled between the boundaries of the
    obj_fun : function, optional
        Function that takes 2 parameters, recorded and simulated discharge. If
        unspecified, RMSE is used.
    wu : int, optional
        Warming up period. This accounts for the number of steps that the model
        is run before calculating the performance function.
    verbose : bool, optional
        If True, displays the result of each model evaluation when performing
        the calibration of the hydrological model.
    tol : float, optional
        Determines the tolerance of the solutions in the optimisaiton process.
    minimise : bool, optional
        If True, the optimisation corresponds to the minimisation of the 
        objective function. If False, the optimial of the objective function is
        maximised.
    fun_nam : str, optional
        Name of the objective function used in calibration. If unspecified, is
        'RMSE'
    
    Returns
    -------
    params : array_like [18]
        Optimal parameter set
    
    performance : float
        Optimal value of the objective function
    '''

    def _cal_fun(par):
        q_sim = simulate(avg_prec[:-1], temp, et, par, p2, init_st=None,
                         ll_temp=None, q_0=10.0)[0]
        if minimise:
            perf = obj_fun(flow[wu:], q_sim[wu:])
        else:
            perf = -obj_fun(flow[wu:], q_sim[wu:])

        if verbose:
            print('{0}: {1}'.format(fun_nam, perf))
        return perf

    # Boundaries
    x_b = zip(x_lb, x_ub)

    # initial guess
    if x_0 is None:
        # Randomly generated
        x_0 = np.random.uniform(x_lb, x_ub)

    # Model optimisation
    par_cal = opt.minimize(_cal_fun, x_0, method='L-BFGS-B', bounds=x_b,
                           tol=tol)
    params = par_cal.x
    performance = par_cal.fun
    return params, performance


def calibrate_sp(flow, avg_prec, temp, et, p2, init_st=None, ll_temp=None,
              x_0=None, x_lb=P_LB, x_ub=P_UB, obj_fun=_rmse, wu=10,
              verbose=False, tol=0.001, minimise=True, fun_nam='RMSE'):
    '''
    
    ## Not tested/implemented yet!
    =========
    Calibrate
    =========

    Run the calibration of the HBV-96. The calibration is used to estimate the
    optimal set of parameters that minimises the difference between 
    observations and modelled discharge.
    
    Parameters
    ----------
    
    flow : array_like [n]
        Measurements of discharge [m3/s]
    avg_prec : array_like [n]
        Average precipitation [mm/h]
    temp : array_like [n]
        Average temperature [C]
    et : array_like [n]
        Potential Evapotranspiration [mm/h] 
    p2 : array_like [2]
        Problem parameter vector setup as:
        [tfac, area]
    init_st : array_like [5], optional
        Initial model states, [sp, sm, uz, lz, wc]. If unspecified, 
        [0.0, 30.0, 30.0, 30.0, 0.0] mm
    ll_temp : array_like [n], optional
        Long term average temptearature. If unspecified, calculated from temp.
    x_0 : array_like [18], optional
        First guess of the parameter vector. If unspecified, a random value
        will be sampled between the boundaries of the 
    x_lb : array_like [18], optional
        Lower boundary of the parameter vector. If unspecified, a random value
        will be sampled between the boundaries of the 
    x_ub : array_like [18], optional
        First guess of the parameter vector. If unspecified, a random value
        will be sampled between the boundaries of the
    obj_fun : function, optional
        Function that takes 2 parameters, recorded and simulated discharge. If
        unspecified, RMSE is used.
    wu : int, optional
        Warming up period. This accounts for the number of steps that the model
        is run before calculating the performance function.
    verbose : bool, optional
        If True, displays the result of each model evaluation when performing
        the calibration of the hydrological model.
    tol : float, optional
        Determines the tolerance of the solutions in the optimisaiton process.
    minimise : bool, optional
        If True, the optimisation corresponds to the minimisation of the 
        objective function. If False, the optimial of the objective function is
        maximised.
    fun_nam : str, optional
        Name of the objective function used in calibration. If unspecified, is
        'RMSE'
    
    Returns
    -------
    params : array_like [18]
        Optimal parameter set
    
    performance : float
        Optimal value of the objective function
    '''

    def _cal_fun(par):
        q_sim = simulate(avg_prec[:-1], temp, et, par, p2, init_st=None,
                         ll_temp=None, q_0=10.0)[0]
        if minimise:
            perf = obj_fun(flow[wu:], q_sim[wu:])
        else:
            perf = -obj_fun(flow[wu:], q_sim[wu:])

        if verbose:
            print('{0}: {1}'.format(fun_nam, perf))
        return perf

    # Boundaries
    x_b = zip(x_lb, x_ub)

    # initial guess
    if x_0 is None:
        # Randomly generated
        x_0 = np.random.uniform(x_lb, x_ub)

    # Model optimisation
    par_cal = opt.minimize(_cal_fun, x_0, method='L-BFGS-B', bounds=x_b,
                           tol=tol)
    params = par_cal.x
    performance = par_cal.fun
    return params, performance

def _sp_test():
#    dem_uri = r'results\DEM_FILL_cut_coarse.tif'  # DEM file
#    dem = gdal.Open(dem_uri)
    flp = gdal.Open(r'results\flow_path_length_coarse.tif')
#    base_dem = gdal.Open(dem_uri)
    shape_base_dem = flp.ReadAsArray().shape
#    sample_dem = base_dem.ReadAsArray()
#    no_val = base_dem.GetRasterBand(1).GetNoDataValue()
#    x_init, xx_span, xy_span, y_init, yy_span, yx_span = base_dem.GetGeoTransform()
    t_steps = 30
    
    sp_prec = np.random.uniform(0, 10, [t_steps, shape_base_dem[0], shape_base_dem[1]])
    sp_temp = np.random.uniform(10, 30, [t_steps, shape_base_dem[0], shape_base_dem[1]])    
    sp_et = np.random.uniform(0.001, 0.01, [t_steps, shape_base_dem[0], shape_base_dem[1]])
    
    sp_pars = P_LB
    p2 = [1, 145.0]
    q_out, st = distributed(flp=flp, sp_prec=sp_prec, sp_et=sp_et, 
                         sp_temp=sp_temp, sp_pars=sp_pars, p2=p2, 
                         init_st=None, ll_temp=None, q_0=None)
#    
if __name__ == '__main__':
    # testing_fun
    import matplotlib.pyplot as plt
#    impo
    
    flp = gdal.Open(r'results\flow_path_length_coarse.tif')
    shape_base_dem = flp.ReadAsArray().shape
    t_steps = 200
    
    sp_prec = np.random.uniform(0, 10, [t_steps, shape_base_dem[0], shape_base_dem[1]])
    sp_temp = np.random.uniform(10, 30, [t_steps, shape_base_dem[0], shape_base_dem[1]])    
    sp_et = np.random.uniform(0.001, 0.01, [t_steps, shape_base_dem[0], shape_base_dem[1]])
    
    
    
    sp_pars = P_LB
    p2 = [1, 145.0]
    q_out, st = distributed(flp=flp, sp_prec=sp_prec, sp_et=sp_et, 
                            sp_temp=sp_temp, sp_pars=sp_pars, p2=p2, 
                            init_st=None, ll_temp=None, q_0=None)
    
