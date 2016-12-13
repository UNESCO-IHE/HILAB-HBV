# HILAB - HBV
This is the first version of the unified HBV model version at UNESCO-IHE. The code is developed under the umbrella of the Hydroinformatics Laboratory at UNESCO-IHE. Currently the model is only implemented in Matlab. Implementation in other languages are foreseen in the near future..

This folder contains the model routines
 - precipitation module <precipitation.m>
 - Snow module <snow.m>
 - soil moisture module <soil.m>
 - response module <response.m>
 - routing module <routing.m>
 
The model is integrated by single time steps (<step_run.m>).
To simulate a complete time series, the <simulate.m> file, can be used.
An example of use is given in <sample_run.m>.

To do:
------

- Code debugging and testing
- Add performance metrics
- Add automatic calibration
- Add sensitivity analysis
- Add data Assimilation
- Add error correctors
