function satellite_scenario_viewer(position_array, t_array, inclination, eccentricity, argument_of_perigee)
    % position_array: 3-by-n matrix of satellite positions in ECI coordinates
    % t_array: array of time values corresponding to the satellite positions
    % inclination: inclination angle in radians
    % eccentricity: eccentricity value
    % argument_of_perigee: argument of perigee in radians

    % Calculate semimajor axis and period
    mu = 398600.436;
    a = nthroot(mu * t_array(end)^2 / (4 * pi^2), 2);
    T = 2 * pi * sqrt(a^3 / mu);
     
    % Calculate mean anomaly
    n = sqrt(mu / a^3);
    M = n * t_array;

    % Calculate eccentric anomaly
    E = kepler_equation(M, eccentricity);

    % Calculate true anomaly
    
    theta = 2 * atan2(sqrt(1 + eccentricity) * sin(E / 2), sqrt(1 - eccentricity) * cos(E / 2));
    M = theta - argument_of_perigee;

    % Calculate satellite scenario
    scenario = struct();
    scenario.latitude = 0;
    scenario.longitude = 0;
    scenario.height = 0;
    scenario.elevationCutoff = 0;
    scenario.day = 1;
    scenario.startTime = 0;
    scenario.periodHours = T / 3600;
    scenario.timeZone = 0;
    scenario.satellite = struct();
    scenario.satellite(1).name = 'Satellite';
    scenario.satellite(1).description = 'Satellite Description';
    scenario.satellite(1).keplerian = struct();
    scenario.satellite(1).keplerian.inclination = rad2deg(inclination);
    scenario.satellite(1).keplerian.eccentricity = eccentricity;
    scenario.satellite(1).keplerian.meanAnomaly = rad2deg(M);
    scenario.satellite(1).keplerian.semimajorAxis = a;
    scenario.satellite(1).keplerian.longitudeAscendingNode = 0;
    scenario.satellite(1).keplerian.argumentOfPerigee = rad2deg(argument_of_perigee);
    scenario.satellite(1).azimuthCutoff = 0;
    scenario.satellite(1).status = 'Operational';

    % Plot satellite scenario viewer
    web(strcat('https://www.gnssplanning.com/#/scenario-viewer/', urlencode(jsonencode(scenario))));
end