# HILAB - HBV
This is the first version of the unified HBV model version at IHE-Delft. The code is developed under the umbrella of the [Hydroinformatics Laboratory](https://hilab.unesco-ihe.org/) at [IHE-Delft](https://www.unesco-ihe.org). Currently the model is implemented in Matlab, Python and Cython. Implementation in other languages are foreseen in the near future.

The Matlab folder contains the model routines
 - precipitation module <precipitation.m>
 - Snow module <snow.m>
 - soil moisture module <soil.m>
 - response module <response.m>
 - routing module <routing.m>
 
 The Python folder contains the model in pure Python and Cython with declaration of static type variables
 - Pure python <HBV96.py>
 - Cython <HBV96x.pyx>
 
The model is integrated by single time steps (<step_run.m> or <HBV96.step_run>).
To simulate a complete time series, the <simulate.m> or <HBV96.simulate>.
An example of use is given in <sample_run.m> for Matlab and <sample_HBV.py> for Python .

To do:
------

- Code debugging and testing
- Add performance metrics (currently only available _nse and _RMSE in python)
- Add automatic calibration (Matlab)
- Add sensitivity analysis
- Add data Assimilation
- Add error correctors
