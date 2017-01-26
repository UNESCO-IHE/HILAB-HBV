function [qNew, spNew, smNew, uzNew, lzNew, wcNew] = step_run(p, p2, v, st)
%{
STEP_RUN This is the main module for the HBV96 (Lindstrom, 1997)

p = Parameter vector
   1. ltt   = Lower temperature limit for snow/rain precipitation [C]
   2. utt   = Upper temperature limit for snow/rain precipitation [C]
   3. ttm   = Temperature limit for melting [C]
   4. cfmax = Degree-day factor
   5. fc    = Field capacity [mm]
   6. ecorr = Evapotranspiration corrector factor
   7. etf   = Temperature corrector factor [C]
   8. lp    = Limit for potential evapotranspiration [mm]
   9. k     = Upper zone recession coefficient
  10. k1    = Lower zone recession coefficient
  11. alpha = Upper zone response coefficient
  12. beta  = Soil moisture parameter
  13. cwh   = Water holding capacity [mm]
  14. cfr   = Refreezing factor
  15. cflux = Maximum capilary flow [mm/hr]
  16. perc  = Percolation [mm/hr]
  17. rfcf  = Rainfall corrector factor
  18. sfcf  = Snowfall corrector factor

p2 = non optimisable parameter vector
  1. tfac = Time factor for unit conversion (Number of hours in time step)
  2. area = Catchment area [Km²]
  
v = inputs
  1. avgPrec = Precipitation [mm]
  2. temp = Temperature [C]
  3. ep = Long terms (monthly) Evapotranspiration [mm]
  4. tm = Long term (monthly) average temperature [C]

st = model states
  1. spOld = Snow pack [mm]
  2. smOld = Soil moisture [mm]
  3. uzOld = Upper zone [mm]
  4. lzOld = Lower zone [mm]
  5. wcOld = Water content in snow pack [mm]

qNew = Outflow [m³/s]
spNew = Snow pack [mm]
smNew = Soil moisture [mm]
uzNew = Upper zone [mm]
lzNew = Lower zone [mm]
wcNew = Water content in snow pack [mm]

%}  

% Parse parameters
  ltt = p(1); % Lower temperature limit for snow/rain precipitation [C]
  utt = p(2); % Upper temperature limit for snow/rain precipitation [C]
  ttm = p(3); % Temperature limit for melting [C]
  cfmax = p(4); % Degree-day factor
  fc = p(5); % Field capacity [mm]
  ecorr = p(6); % Evapotranspiration corrector factor
  etf = p(7); % Temperature corrector factor [C]
  lp = p(8); % Limit for potential evapotranspiration [mm]
  k = p(9); % Upper zone recession coefficient
  k1 = p(10); % Lower zone recession coefficient
  alpha = p(11); % Upper zone response coefficient
  beta = p(12); % Soil moisture parameter
  cwh = p(13); % Water holding capacity [mm]
  cfr = p(14); % Refreezing factor
  cflux = p(15); % Maximum capilary flow [mm/hr]
  perc = p(16); % Percolation [mm/hr]
  rfcf = p(17); % Rainfall corrector factor
  sfcf = p(18); % Snowfall corrector factor

  tfac = p2(1); % Time factor for unit conversion (Number of hours in time step)
  area = p2(2); % Catchment Area [Km²]
  
  % selecting the snow switch
  if length(p2) == 3;
    snow_switch = p2(3);
  else
    snow_switch = 1;
  endif

% Parse inputs  
avgPrec = v(1); % Precipitation [mm]
temp = v(2); % Temperature [C]
ep = v(3); % Long terms (monthly) Evapotranspiration [mm]
tm = v(4); %Long term (monthly) average temperature [C]

% Parse states
spOld = st(1); % Snow pack [mm]
smOld = st(2); % Soil moisture [mm]
uzOld = st(3); % Upper zone [mm]
lzOld = st(4); % Lower zone [mm]
wcOld = st(5); % Water content in snow pack [mm]

if snow_switch == 1;
  % Precipitation routine
  [rf, sf] = precipitation(temp, ltt, utt, avgPrec, rfcf, sfcf, tfac);
  
  % Snow routine
  [inf, wcNew, spNew]	= snow(cfmax, tfac, temp, ttm, cfr, cwh, rf, ...
                                sf, wcOld, spOld);
else
  inf = avgPrec;
  wcNew = 0;
  spNew = 0;

endif

% Soil moisture routine 
[smNew, uzInt_1, qDr] = soil(fc, beta, etf, temp, tm, ecorr, lp,...
                               tfac, cflux, inf, ep, smOld, uzOld);

% Response routine
[qNew, uzNew, lzNew] = response(tfac, perc, alpha, k, k1, area, ...
                                   lzOld, uzInt_1, qDr);

end
