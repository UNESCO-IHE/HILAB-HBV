
function [smNew, uzInt_1, qDr] = soil(fc, beta, etf, temp, tm, ecorr,...
                                        lp, tfac, cflux, inf, et, ...
                                        smOld, uzOld)
  %{
  SOIL runs the soild moisture routine of the HBV96

  At first, comparison of temperature is made. if temperature is below
  treshold, melting is happening, otherwise, refreezing. if the water
  content in the snow pack is bigger than water holding capacity, excess
  infiltrates soil.

       fc = Filed capacity [mm]
     beta = Shape coefficient for effective precipitation separation
      etf = Total potential evapotranspiration [mm]
      temp = Temperature [C]
       tm = Average long term temperature [C]
    ecorr = Evapotranspiration corrector factor
       lp = Evapotranspiration limit factor
     tfac = Time conversion factor [mm]
    cflux = Capilar flux limit in the root zone [mm/hr]
      inf = actual infiltration [mm]
       et = evapotranspiration [mm]
    smOld = Previous soil moisture value [mm]
    uzOld = Previous Upper zone value [mm]
  
    smNew = New value of soil moisture [mm]
  uzInt_1 = New value of direct runoff into upper zone [mm]
      qDr = Direct runoff [mm]
    %}

    qDr = max(smOld + inf - fc, 0); % Overland discharge
    in = inf - qDr; % Effective infiltration
    r = ((smOld/fc)^beta) * in; % Direct runoff
    epInt = max((1.0 + etf*(temp - tm))* ecorr*et, 0); % day degree
    ea = max(epInt, (smOld/(lp*fc))* epInt); % Actual evapotranspiration
    cf = min(cflux*((fc - smOld)/fc), uzOld + r); % Capilar flow
    smNew = max(smOld + in - r + cf - ea, 0); # Soil moisture update
    uzInt_1 = uzOld + r - cf; % percolation to UZ
end
