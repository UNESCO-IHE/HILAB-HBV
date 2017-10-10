function [qNew, uzNew, lzNew] = response(tFac, perc, alpha, k, k1, ...
                                            area, lzOld, uzInt_1, qDr)
  % RESPONSE makes the catchment response (upper and lower zone)of the HBV96 
  % model
  % 
  % [q_new, uz_new, lz_new] = response(tfac, perc, alpha, k, k1, ...
  %                                    area, lz_old, uz_int_1, qdr)
  % 
  % This solution has preferential recharge over overland flow
  %     tFac = number of hours in the time step
  %     perc = Percolation from upper to lower response box [mm/hr]
  %    alpha = Response box parameter
  %        k = Upper zone recession coefficient
  %       k1 = Lower zone recession coefficient
  %     area = Catchment area [km²]
  %   lzOld = Lower zone
  % uzInt_1 = Upper zone after soil moisture update
  %      qDr = Direct runoff [mm]
  %
  %    qNew = discharge before routing [m³/s]
  %   uzNew = updated upper zone [mm]
  %   lzNew = updated lower zone [mm]

  % calculate changes in the lower and upper zones by percolation
  lzInt_1 = lzOld + min(perc, uzInt_1);
  uzInt_2 = max(uzInt_1 - perc, 0);

  % Calculate runoff from upper, lower zone and direct runoff
  q0 = min(k*(uzInt_2^(1.0 + alpha)), uzInt_2);
  q1 = min(k1*lzInt_1, lzInt_1);
  qNew = area*(q0 + q1 + qDr)/(tFac*3.6);

  % Update upper and lower zone
  uzNew = max(uzInt_2 - (q0), 0);
  lzNew = max(lzInt_1 - (q1), 0);


end
