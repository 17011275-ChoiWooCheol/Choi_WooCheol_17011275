function az = azimuth(ENU)
    % ENU: n-by-3 matrix of satellite ENU positions (km)
    
    % Calculate RN (Satellite's distance from the observer in the North direction)
    RN = ENU(:, 2);
    
    % Calculate RE (Satellite's distance from the observer in the East direction)
    RE = ENU(:, 1);
    
    % Calculate azimuth angle
    az = acosd(RN ./ sqrt(RE.^2 + RN.^2));
end