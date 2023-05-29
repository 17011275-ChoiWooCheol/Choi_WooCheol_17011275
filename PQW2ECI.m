function R_pqw_eci = PQW2ECI(arg_prg, inc_angle, RAAN)
    %-degree to radian
    arg_prg = deg2rad(arg_prg);
    inc_angle = deg2rad(inc_angle);
    RAAN = deg2rad(RAAN);

    % rotation matrix elements
    cos_arg_prg = cos(arg_prg);
    sin_arg_prg = sin(arg_prg);
    cos_inc_angle = cos(inc_angle);
    sin_inc_angle = sin(inc_angle);
    cos_RAAN = cos(RAAN);
    sin_RAAN = sin(RAAN);

    % Construct the rotation matrix
    R_pqw_eci = [-sin_RAAN*cos_inc_angle*sin_arg_prg + cos_RAAN*cos_arg_prg, ...
                 -sin_RAAN*cos_inc_angle*cos_arg_prg - cos_RAAN*sin_arg_prg, ...
                 sin_RAAN*sin_inc_angle;
                 cos_RAAN*cos_inc_angle*sin_arg_prg + sin_RAAN*cos_arg_prg , ...
                 cos_RAAN*cos_inc_angle*cos_arg_prg - sin_RAAN*sin_arg_prg , ...
                 -cos_RAAN*sin_inc_angle;   
                 sin_inc_angle*sin_arg_prg sin_inc_angle*cos_arg_prg cos_inc_angle];
end
