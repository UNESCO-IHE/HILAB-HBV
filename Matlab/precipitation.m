function [rf, sf] = precipitation(temp, ltt, utt, prec, rfcf, sfcf, tfac)
    %{
    Precipitation routine
    ---------------------
    
    Below the lower treshold level, all precipitation is snowfall.
    Similarly, all the precipitation above upper temperature treshold is
    rainfall. In between, there is a linear mixture between raifall and
    snowfall.
    
    Parameters
    ----------
        **temp -- Measured temperature (C)
        **ltt -- Lower temperature treshold (C)
        **utt -- Upper temperature treshold (C)
        **prec -- Precipitation (mm)
        **rfcf -- Rainfall corrector factor
        **sfcf -- Snowfall corrector factor
    Returns
    -------
        **_rf - Rainfall (mm)
        **_sf - Snowfall (mm)
    %}
    if temp <= ltt;
        rf = 0.0;
        sf = prec*sfcf;
		
    elseif temp >= utt;
        rf = prec*rfcf;
        sf = 0.0;
        
    else;
        rf = ((temp-ltt)/(utt-ltt)) * prec * rfcf;
        sf = (1.0-((temp-ltt)/(utt-ltt))) * prec * sfcf;
    end
end
