function wx1()
    % 주어진 궤도 정보
    semi_major_axis = 6400;                         % 반지름
    eccentricity = 0.0829;                          % 이심률
    inclination = deg2rad(30);                      % 경사각
    argument_of_perigee = deg2rad(60);              % 근지점의 위상
    right_ascension = deg2rad(45);                  % 승교점의 적경
    mean_anomaly = deg2rad(0);                      % 평균 근점 이각
    t_start = 0;                                    % 시작 시각
    t_end = 24 * 3600;                              % 종료 시각

    % 중력 상수
    mu = 398600.436;

    % 시간 단계
    dt = 60;

    % 시간 배열 생성
    t_array = t_start:dt:t_end;

    % 위치 배열 초기화
    position_array = zeros(3, numel(t_array));

    % 시간에 따른 위성 위치 계산
    for i = 1:numel(t_array)
        t = t_array(i);
        [position, ~] = keplerian_to_cartesian(semi_major_axis, eccentricity, inclination, argument_of_perigee, right_ascension, mean_anomaly, t);
        position_array(:, i) = position;
    end

    % 위성 궤도 그리기
    figure;
    plot3(position_array(1, :), position_array(2, :), position_array(3, :));
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Satellite Orbit');

    % 지상 궤적 그리기
    figure;
    ground_track = sqrt(position_array(1, :).^2 + position_array(2, :).^2);
    plot(t_array, ground_track);
    xlabel('Time (s)');
    ylabel('Ground Track');
    title('Satellite Ground Track');
    %%position_array = zeros(3, numel(t_array));
    % skyplot 그리기
    figure;
    azi = atan2(position_array(2, :), position_array(1, :));
    ele = asin(position_array(3, :) ./ sqrt(sum(position_array.^2)));
    polarplot(azi, 90 - rad2deg(ele), '.');
    title('Satellite Skyplot');

    % sky view 그리기
        figure;
        sky_view = zeros(181, 361);
        for i = 1:numel(azi)
            az = rad2deg(azi(i));
            el = rad2deg(ele(i));
            az_idx = round(az) + 181; % 인덱스 수정
            el_idx = round(el) + 91; % 인덱스 수정
            sky_view(el_idx, az_idx) = 1;
        end
        [X, Y] = meshgrid(-180:180, -90:90);
        pcolor(X, Y, sky_view);
        shading flat;
        colormap(gray);
        colorbar;
        caxis([0 1]);
        xlabel('Azimuth (degrees)');
        ylabel('Elevation (degrees)');
        title('Sky View');

    % Calculate true anomaly
    E = 2 * atan(sqrt((1 - eccentricity) / (1 + eccentricity)) * tan(mean_anomaly / 2));

    % satellite scenario viewer 호출
    satellite_scenario_viewer(position_array, t_array, inclination, eccentricity, argument_of_perigee);
end

function [position, velocity] = keplerian_to_cartesian(semi_major_axis, eccentricity, inclination, argument_of_perigee, right_ascension, mean_anomaly, t)
    % 중력 상수
    mu = 398600.436;

    % 평균 운동 계산
    n = sqrt(mu / semi_major_axis^3);

    % 평균 운동 anomaly 계산
    M = mean_anomaly + n * t;

    % 반복적인 평균 운동 anomaly 계산
    while M < 0
        M = M + 2 * pi;
    end

    % 이심 운동 anomaly 계산
    E = kepler_equation(M, eccentricity);

    % 진근점 거리 계산
    r = semi_major_axis * (1 - eccentricity * cos(E));

    % 위경도 계산
    nu = atan2(sqrt(1 - eccentricity^2) * sin(E), cos(E) - eccentricity);

    % 위성 위치 계산
    position = [
        r * (cos(right_ascension) * cos(argument_of_perigee + nu) - sin(right_ascension) * sin(argument_of_perigee + nu) * cos(inclination));
        r * (sin(right_ascension) * cos(argument_of_perigee + nu) + cos(right_ascension) * sin(argument_of_perigee + nu) * cos(inclination));
        r * sin(argument_of_perigee + nu) * sin(inclination)
    ];

    % 위성 속도 계산
    velocity = [
        -mu * sin(argument_of_perigee + nu) / r;
        mu * (eccentricity + cos(argument_of_perigee + nu)) / r;
        0
    ];
end

function E = kepler_equation(M, eccentricity)
    E = M;                                       % 초기값 설정
    threshold = 1e-8;                            % 수렴 조건 임계값
    error = 1;                                   % 오차 초기값 설정

    % 반복적으로 이심 운동 방정식 해 구하기
    while error > threshold
        E_new = E - (E - eccentricity * sin(E) - M) / (1 - eccentricity ...
            * cos(E));
        error = abs(E_new - E);
        E = E_new;
    end
end

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
    nu = 2 * atan(sqrt((1 + eccentricity) / (1 - eccentricity)) * tan(E / 2));

    % Calculate argument of latitude
    u = nu + argument_of_perigee;

    % Calculate longitude of ascending node
    Omega = 0;

    % Calculate satellite positions in ECI coordinates
    X = a * (cos(u) * cos(Omega) - sin(u) * sin(Omega) * cos(inclination));
    Y = a * (cos(u) * sin(Omega) + sin(u) * cos(Omega) * cos(inclination));
    Z = a * (sin(u) * sin(inclination));

    % Plot satellite positions
    figure;
    plot3(X, Y, Z);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Satellite Scenario Viewer');
end
    