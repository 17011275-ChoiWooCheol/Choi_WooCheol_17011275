function el = elevation(ENU, el_mask)
    % ENU: n-by-3 matrix of satellite ENU positions in km
    % el_mask: Elevation mask angle (deg)
    
    % Calculate Satellite's distance from the observer in the up direction
    Ru = ENU(:, 3);
    
    % Calculate Relative distance from the observer to the satellite
    Rrel = vecnorm(ENU, 2, 2);
    
    % Calculate elevation angle
    el = asind(Ru ./ Rrel);
    
    el(el < el_mask) = NaN;
end