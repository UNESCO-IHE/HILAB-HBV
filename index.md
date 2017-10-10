# HILAB - HBV
This is the first version of the unified HBV model version at UNESCO-IHE. The code is developed under the umbrella of the [Hydroinformatics Laboratory](https://hilab.unesco-ihe.org/) at [UNESCO-IHE](https://www.unesco-ihe.org). Currently the model is only implemented in Matlab. Implementation in other languages are foreseen in the near future.

This folder contains the model routines
 - precipitation module (precipitation)
 - Snow module (snow)
 - soil moisture module (soil)
 - response module (response)
 - routing module (routing)
 
The model is integrated by single time steps (`step_run.m`).
To simulate a complete time series, he `simulate.m` file, can be used.
A running example of use is provided in `sample_run.m`

## To do:

- Code debugging and testing
- Add performance metrics
- Add automatic calibration
- Add sensitivity analysis
- Add data assimilation
- Add error corrector
