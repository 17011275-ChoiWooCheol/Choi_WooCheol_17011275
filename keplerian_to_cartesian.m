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
    E = M; % 초기값 설정
    threshold = 1e-8; % 수렴 조건 임계값
    error = 1; % 오차 초기값 설정

    % 반복적으로 이심 운동 방정식 해 구하기
    while error > threshold
        E_new = E - (E - eccentricity * sin(E) - M) / (1 - eccentricity * cos(E));
        error = abs(E_new - E);
        E = E_new;
    end
end
