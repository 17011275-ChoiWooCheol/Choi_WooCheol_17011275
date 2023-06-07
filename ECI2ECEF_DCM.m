function DCM = ECI2ECEF_DCM(time)
    % time: [YYYY, MM, DD, hh, mm, ss] 
    % Convert time to serial date number
    t = datetime(time);
    serial_date = datenum(t);
    
    % Calculate Julian Date
    JD = serial_date - datenum([2000, 1, 1, 12, 0, 0]);
    
    % Calculate GMST in radians
    GMST = deg2rad(GMST_JD(JD));
    
    % Calculate rotation angle for ECI to ECEF
    theta = GMST + deg2rad(360) * (JD - 1);
    
    % Construct DCM Matrix
    DCM = [cos(theta), sin(theta), 0;
           -sin(theta), cos(theta), 0;
           0, 0, 1];
end