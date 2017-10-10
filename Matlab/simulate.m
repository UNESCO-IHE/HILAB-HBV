function [qRout, st] = simulate(avgPrec, temp, et, par, p2)
%{
SIMULATE runs the HBV96 for the given dataset

avgPrec = Precipitation [mm]
temp = Temperature [C]
et = Long terms (monthly) Evapotranspiration [mm]

par = Parameter vector
   1. ltt    = Lower temperature limit for snow/rain precipitation [C]
   2. utt    = Upper temperature limit for snow/rain precipitation [C]
   3. ttm    = Temperature limit for melting [C]
   4. cfmax  = Degree-day factor
   5. fc     = Field capacity [mm]
   6. ecorr  = Evapotranspiration corrector factor
   7. etf    = Temperature corrector factor [C]
   8. lp     = Limit for potential evapotranspiration [mm]
   9. k      = Upper zone recession coefficient
  10. k1     = Lower zone recession coefficient
  11. alpha  = Upper zone response coefficient
  12. beta   = Soil moisture parameter
  13. cwh    = Water holding capacity [mm]
  14. cfr    = Refreezing factor
  15. cflux  = Maximum capilary flow [mm/hr]
  16. perc   = Percolation [mm/hr]
  17. rfcf   = Rainfall corrector factor
  18. sfcf   = Snowfall corrector factor
  19. maxbas = Flow routing coefficient

p2 = non optimisable parameter vector
  1. tfac = Time factor for unit conversion (Number of hours in time step)
  2. area = Catchment area [Km²]

%}

  % Initialize model run
  q0 = 10.0;
  st = zeros(length(avgPrec), 5);
  st(1,:) = [30, 30, 30, 30, 30];
  llTemp = mean(temp) * ones(length(avgPrec), 1);
  qSim = zeros(length(avgPrec), 1);
  qSim(1) = q0;

  % Run the model for the available data
  for i=1:length(avgPrec)-1;
    v = [avgPrec(i), temp(i), et(i), llTemp(i)];
    [qOut, stOut(1), stOut(2), stOut(3), stOut(4), stOut(5)] = step_run(par, p2, v, st(i,:));
    qSim(i+1) = qOut;
    st(i+1,:) = [stOut(1), stOut(2), stOut(3), stOut(4), stOut(5)];
  end

  % Do the routing if maxbas is included
  qRout = qSim;
  if par(end) >= 2;
    qRout = routing(qSim, par(end));  
  end

end


