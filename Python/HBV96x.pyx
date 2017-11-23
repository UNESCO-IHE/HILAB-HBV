# -*- coding: utf-8 -*-
"""
Created on Mon Jun 03 14:17:49 2013

@author: chaco3

HBV-96
This is the HBV-96 implementation by Juan Chacon at UNESCO-IHE, Delft, NL.
This is meant to be part of the wesenseit project,of citizens observatories
"""
from __future__ import division, print_function
import numpy as np
import scipy.optimize as opt

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

# Get random parameter set
def get_random_pars():
    return np.random.uniform(P_LB, P_UB)

def _precipitation(double temp, double ltt, double utt, double prec, 
                   double rfcf, double sfcf, double tfac):
    '''
    Precipitation routine
    ---------------------
    
    Below the lower treshold level, all precipitation is snowfall.
    Similarly, all the precipitation above upper temperature treshold is
    rainfall. In between, there is a linear mixture between raifall and
    snowfall.
    
    Parameters
    ----------
        **temp -- Measured temperature [C]
        **ltt -- Lower temperature treshold [C]
        **utt -- Upper temperature treshold [C]
        **prec -- Precipitation [mm]
        **rfcf -- Rainfall corrector factor
        **sfcf -- Snowfall corrector factor
    Returns
    -------
        **_rf - Rainfall [mm]
        **_sf - Snowfall [mm]
    '''
    
    cdef double _rf,_sf
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


def _snow(double cfmax, double tfac, double temp, double ttm, double cfr, 
          double cwh, double _rf, double _sf, double wc_old, double sp_old):
    '''
    _snow routine
    ------------
        At first, comparison of temperature is made. if temperature is below
        treshold, melting is happening, otherwise, refreezing. if the water
        content in the snow pack is bigger than water holding capacity, excess
        infiltrates soil.
        
    Parameters
    ----------
        **cfmax -- Day degree factor
        **tfac -- Temperature correction factor
        **temp -- Temperature [C]
        **ttm -- Temperature treshold for Melting [C]
        **cfr -- Refreezing factor
        **cwh -- Capacity for water holding in snow pack
        **_rf -- Rainfall [mm]
        **_sf -- Snowfall [mm]
        **wc_old -- Water content in previous state [mm]
        **sp_old -- _snow pack in previous state [mm]
        
    Returns
    -------
        **_in -- Infiltration [mm]
        **_wc_new -- Water content in posterior state [mm]
        **_sp_new -- Snowpack in posterior state [mm]   
    '''
             
    cdef double _melt, _sp_new, _wc_int, _refr, _in, _wc_new
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


def _soil(double fc, double beta, double etf, double temp, double tm,
          double e_corr, double lp, double tfac, double c_flux, double inf, 
          double ep, double sm_old, double uz_old):
    '''
    _soil routine
    ------------
        At first, comparison of temperature is made. if temperature is below
        treshold, melting is happening, otherwise, refreezing. if the water
        content in the snow pack is bigger than water holding capacity, excess
        infiltrates soil.
        
    Parameters
    ----------
        **fc -- Filed capacity
        **beta -- Shape coefficient for effective precipitation separation
        **etf -- Total potential evapotranspiration
        **temp -- Temperature
        **tm -- Average long term temperature
        **e_corr -- Evapotranspiration corrector factor
        **lp -- _soil wilting point
        **tfac -- Time conversion factor
        **c_flux -- Capilar flux in the root zone
        **_in -- actual infiltration
        **ep -- actual evapotranspiration
        **sm_old -- Previous soil moisture value
        **uz_old -- Previous Upper zone value
        
    Returns
    -------
        **sm_new -- New value of soil moisture
        **uz_int_1 -- New value of direct runoff into upper zone
    '''
    
    cdef double _r, _ep_int, _ea, _cf, sm_new, uz_int_1, qdr
#    ep = ep/tfac
    qdr = max(sm_old + inf - fc, 0)
    _in = inf - qdr
    _r = ((sm_old/fc)** beta) * _in
    _ep_int = (1.0 + etf*(temp - tm))*e_corr*ep
    
#    _r = ((sm_old/fc) ** beta) * _in
#    _ep_int = e_corr*ep
#    
#    if sm_old/(lp*fc) < 1.0: 
#        _ea = (sm_old/(lp*fc))*_ep_int 
#    else: 
#        _ea = _ep_int 
        
    _ea = max(_ep_int, (sm_old/(lp*fc))*_ep_int)
#    if c_flux*(1.0 - (sm_old/fc)) < uz_old:
#        _cf = c_flux*(1.0 - (sm_old/fc))
#    else:
#        _cf = uz_old 
    _cf = c_flux*((fc - sm_old)/fc)
    sm_new = max(sm_old + _in - _r + _cf - _ea, 0)        
    uz_int_1 = uz_old + _r - _cf
    
    return sm_new, uz_int_1, qdr


def _response(double tfac, double perc, double alpha, double k, double k1, 
              double area, double lz_old, double uz_int_1, double qdr):
             
    cdef double lz_int_1, uz_int_2, _q_0, _q_1, uz_new, lz_new, q_new
    
    if perc < uz_int_1: 
        lz_int_1 = lz_old + perc
    else:
        lz_int_1 = lz_old + uz_int_1 

    if uz_int_1 > perc: 
        uz_int_2 = uz_int_1 - perc
    else:
        uz_int_2 = 0.0

    _q_0 = k*(uz_int_2**(1.0 + alpha))
    _q_1 = k1*lz_int_1
    
    uz_new = max(uz_int_2 - (_q_0), 0)
    lz_new = max(lz_int_1 - (_q_1), 0)
    
    q_new = area*(_q_0 + _q_1 + qdr)/(3.6)
    return q_new, uz_new, lz_new


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
        q_r += np.array(q_temp)*w_i
        q_temp = np.insert(q_temp, 0, 0.0)[:-1]

    return q_r
"""
## Routing routine

#  Data load
#  v = Input vector
#  p = Parameter vector
#  St = State vector
#  x = Aditional parameters (values that multiply old states in the model)

# Input variables 
# prec = Total precipitation v(0)
# temp = actual temperature v(1)
# etf = Total potential evapotranspiration v(2)   input
# tm = daily long term mean temperature v(3)  input

# Parameter Set
# TT = Limit temperature for rain/snow precipitation p(0)
# TTI = temperature treshold for linear mix of snow/rain precipitation p(1)
# ttm = Limit temperature for melting p(2)
# cfmax = Degree day factor (measures the temperature variation along the day) p(3)
# fc = Field Capacity p(4)
# e_corr = Evapotranspiration corrector factor p(5)
# ep = Long term mean potential evapotranspiration p(6)
# lp = _soil moisture value where soil moisture reaches maximum potential
# evapotranspiration p(7)
# k = Upper zone response coefficient p(8)
# k1 = Lowe zone response coefficient p(9)
# alpha = upper zone runoff coefficient p(10)
# beta = Controls the contribution of the increase in the soil moisture or
# to the response function p(11)
# cwh = Maximum amount of water that can be stored in snow pack p(12)
# cfr = Refreezing factor p(13)
# c_flux = Capilar flux p(14)
# perc = Percolation p(15)
# rfcf = Rainfal correction factor p(16)
# sfcf = Snowfall correction factor p(17)

# Non optimised parameters
# tfac = Time factor p(18) = dt/86400
# area = Catchment area p(19)
"""


def _step_run(p, p2, v, St):
    '''
    This is the main module for the HBV run, as described in Lindstrom, 1997
    
    This script receives
        p = Parameter vector
        p2 = non optimisable parameter vector
        v = inputs
        St = Old states of the model
    
    This script returns
        q_new = Outflow
        St = Posterior states of the model
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
    sm_new, uz_int_1, qdr = _soil(fc, beta, etf, temp, tm, e_corr, lp, tfac, c_flux, 
                            inf, ep, sm_old, uz_old)
    q_new, uz_new, lz_new = _response(tfac, perc, alpha, k, k1, area, lz_old, 
                                    uz_int_1, qdr)
    
    return q_new, [sp_new, sm_new, uz_new, lz_new, wc_new]


def simulate(avg_prec, temp, et, par, p2, init_st=None, ll_temp=None, q_0=10.0):
    '''
    ========
    Simulate
    ========
    
    Run the HBV model for the number of steps (n) in precipitation. The
    resluts are (n+1) simulation of discharge as the model calculates step n+1
    
    
    '''    
    if init_st is None:
        st = [[30.0, ]*5]
    else:
        st = [init_st, ]

    if ll_temp is None:
        ll_temp = [np.mean(temp), ] * len(avg_prec)
    
    q_sim = [q_0, ]

    for i in xrange(len(avg_prec)):
        v = [avg_prec[i], temp[i], et[i], ll_temp[i]]
        q_out, st_out = _step_run(par, p2, v, st[i])
        q_sim.append(q_out)
        st.append(st_out)
    
    if len(p2) > 2:  # Forcing maxbas to be predefined
        maxbas = p2[2]  
    elif len(par) > 18:  # Putting maxbas as parameter to be optimised
        maxbas = par[18]
    else:
        maxbas = 1
        
    q_tr = _routing(q_sim, maxbas)
    
    return q_tr, st

def _nse(x, y):
    '''
    x = measured
    y - Simulated
    '''
    a = np.square(np.subtract(x, y))
    b = np.square(np.subtract(x, np.nanmean(x)))
    if a.any < 0.0:
        return(np.nan)
    cdef double f = 1.0 - (np.nansum(a)/np.nansum(b))
    return(f)


def _rmse(x,y):
    '''
    x = measured
    y - Simulated
    '''
    erro = np.square(np.subtract(x,y))
    if erro.any < 0:
        return(np.nan)
    cdef double f = np.sqrt(1.*np.nanmean(erro))
    return(f)

def calibrate(flow, avg_prec, temp, et, p2, init_st=None, ll_temp=None,
              x_0=None, x_lb=P_LB, x_ub=P_UB, obj_fun=_rmse, wu=10,
              verbose=False, tol=0.001, minimise=True, fun_nam='RMSE'):
    '''
    =========
    Calibrate
    =========
    
    Running the calibration of the HBV-96
    '''
    
    if obj_fun == 'NSE':
        def _cal_fun(par):
            q_sim = simulate(avg_prec[:-1], temp, et, par, p2, init_st=None, 
                             ll_temp=None, q_0=10.0)[0]
            perf = -_nse(flow[wu:], q_sim[wu:])
            if verbose:
                print('NSE: {0}'.format(-perf))
            return perf
            
    if obj_fun == 'RMSE':
        def _cal_fun(par):
            q_sim = simulate(avg_prec[:-1], temp, et, par, p2, init_st=None, 
                             ll_temp=None, q_0=10.0)[0]
            perf = _rmse(flow[wu:], q_sim[wu:])
            if verbose:
                print('RMSE: {0}'.format(perf))
            return perf
        
    # initial guess
    if x_0 is None:
        # Randomly generated
        x_0 = np.random.uniform(x_lb, x_ub)
    
    # Boundaries
    if x_lb is None:
        x_lb = P_LB
    
    if x_ub is None:
        x_ub = P_UB
        
    x_b = zip(x_lb, x_ub)
    
    # Model optimisation
    par_cal = opt.minimize(_cal_fun, x_0, method='L-BFGS-B', bounds=x_b, tol=tol)
    params = par_cal.x
    performance = par_cal.fun
    return params, performance


