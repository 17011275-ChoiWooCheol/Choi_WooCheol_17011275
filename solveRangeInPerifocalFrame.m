function rangeInPQW = solveRangeInPerifocalFrame(semimajor_axis,eccentricity,true_anomaly)
       
       p = semimajor_axis*(1 - eccentricity^2); %Semi-Latus-Rectum

       r = p / (1 + (eccentricity * cosd(true_anomaly))); %range in scalar
       rangeInPQW = [r * cosd(true_anomaly); r * sind(true_anomaly); 0]; %3 by 1 matirx
end


