function velocityInPQW = solveVelocityInPerifocalFrame(semimajor_axis,eccentricity, true_anomaly)
        
        mu = 3.986004418*10^14; % Gravitational Parameter 
        p = semimajor_axis*(1 - eccentricity^2); %Semi-Latus-Rectum
        
        A = sqrt(mu/p); %compute in one-parameter
        
        velocityInPQW = [A*(-sind(true_anomaly)); A*(eccentricity + cosd(true_anomaly)); 0];
end