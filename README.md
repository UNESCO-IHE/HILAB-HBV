# MLAB-HBV96
This is the first version of the unified HBV96 version at UNESCO-IHE. This first version is based in the original code by Juan Chacon (https://www.unesco-ihe.org/juan-carlos-chacon-hurtado).

This folder contains the model routines
 - precipitation module (precipitation)
 - Snow module (snow)
 - soil moisture module (soil)
 - response module (response)
 - routing module (routing)
 
The model is integrated by single time steps (step_run.m).
To simulate a complete time series, the simulate.m file, can be used.
An example of use is given in sample_run.m

To do:
------

- Code debugging and testing
- Add performance metrics
- Add automatic calibration
- Add sensitivity analysis
- Add data Assimilation
- Add error correctors
