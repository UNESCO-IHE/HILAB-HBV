function [in, wc_new, sp_new] = snow(cfmax, tfac, temp, ttm, cfr, cwh, rf, sf, wc_old, sp_old)
    %{
    snow routine
    ------------
        At first, comparison of temperature is made. if temperature is below
        treshold, melting is happening, otherwise, refreezing. if the water
        content in the snow pack is bigger than water holding capacity, excess
        infiltrates soil.
        
    Parameters
    ----------
        **cfmax -- Day degree factor
        **tfac -- Temperature correction factor
        **temp -- Temperature (C)
        **ttm -- Temperature treshold for Melting (C)
        **cfr -- Refreezing factor
        **cwh -- Capacity for water holding in snow pack
        **_rf -- Rainfall (mm)
        **_sf -- Snowfall (mm)
        **wc_old -- Water content in previous state (mm)
        **sp_old -- snow pack in previous state (mm)
        
    Returns
    -------
        **_in -- Infiltration (mm)
        **_wc_new -- Water content in posterior state (mm)
        **_sp_new -- Snowpack in posterior state (mm)   
    %}
       
    if temp > ttm;
        if cfmax*(temp-ttm) < sp_old+sf;
            melt = cfmax*(temp-ttm);
        else
            melt = sp_old + sf;
        end
        sp_new = sp_old + sf - melt;
        wc_int = wc_old + melt + rf;

    else
        if cfr*cfmax*(ttm-temp) < wc_old+rf;
            refr = cfr*cfmax*(ttm - temp);
        else
            refr = wc_old + rf;
        end 
        sp_new = sp_old + sf + refr;
        wc_int = wc_old - refr + rf;
		end
	
    if wc_int > cwh*sp_new;
        in = wc_int-cwh*sp_new;
        wc_new = cwh*sp_new;
    else
        in = 0.0;
        wc_new = wc_int;
    end
end
