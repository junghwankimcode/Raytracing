clear; clc; close all;

fc = 2.4e9; 
num_rx = 1e3; 
min_radius = 10; % (m)
max_radius = 50; % (m)
tx1_power_dbm = 10; % Tx1 power (dBm)
tx2_power_dbm = 30; % Tx2 power (dBm)
noise_psd_dbm_per_hz = -174; % dBm/Hz
path_loss_exponent = 3; 
B = 20e6; % Bandwidth (Hz)
noise_variance_watt = 10^((noise_psd_dbm_per_hz - 30) / 10) * B; % Noise power in watts

epsilon_r_real_vec = [5, 10, 15]; % Relative permittivity for metasurfaces
conductivity_vec = [1e-1]; % Conductivity for typical RIS materials

sinr_results = zeros(num_rx, length(epsilon_r_real_vec), length(conductivity_vec));

% Set up Tx sites
tx1 = txsite(Name="Tx1", Latitude=37.5584876, Longitude=126.9997299, ...
    AntennaHeight=3, TransmitterPower=tx1_power_dbm, TransmitterFrequency=fc);

tx2 = txsite(Name="Tx2", Latitude=37.5582876, Longitude=127.0009671, ...
    AntennaHeight=3, TransmitterPower=tx2_power_dbm, TransmitterFrequency=fc);

rx_positions = cell(1, num_rx);
max_cir_length = 0;

% Generate Rx positions
for i = 1:num_rx
    r = sqrt(rand() * (max_radius^2 - min_radius^2) + min_radius^2);
    theta = 2 * pi * rand();
    rx_latitude = tx1.Latitude + (r / 111000) * cos(theta);
    rx_longitude = tx1.Longitude + (r / (111000 * cosd(tx1.Latitude))) * sin(theta);
    rx_positions{i} = rxsite(Name="RX", Latitude=rx_latitude, Longitude=rx_longitude, AntennaHeight=1.7);

    pm1 = propagationModel("raytracing", MaxNumReflections=2, AngularSeparation="medium", ...
        Method="image", BuildingsMaterial="perfect-reflector", SurfaceMaterial="plasterboard");
    pm2 = propagationModel("raytracing", MaxNumReflections=2, AngularSeparation="medium", ...
        Method="image", BuildingsMaterial="concrete", SurfaceMaterial="plasterboard");

    rays_tx1 = raytrace(tx1, rx_positions{i}, pm1);
    rays_tx2 = raytrace(tx2, rx_positions{i}, pm2);

    if ~isempty(rays_tx1)
        pathDelays_tx1 = [rays_tx1{1, 1}.PropagationDelay];
        max_cir_length = max(max_cir_length, ceil(max(pathDelays_tx1) * 1e9) + 1);
    end
    if ~isempty(rays_tx2)
        pathDelays_tx2 = [rays_tx2{1, 1}.PropagationDelay];
        max_cir_length = max(max_cir_length, ceil(max(pathDelays_tx2) * 1e9) + 1);
    end
end

t = 0:1e-9:(max_cir_length - 1) * 1e-9;

% Main simulation loop
epsilon_0 = 8.854e-12; % Vacuum permittivity (F/m)
for e_idx = 1:length(epsilon_r_real_vec)
    for c_idx = 1:length(conductivity_vec)
        epsilon_r_real = epsilon_r_real_vec(e_idx);
        conductivity = conductivity_vec(c_idx);
        complex_permittivity = epsilon_r_real - 1j * conductivity / (2 * pi * fc * epsilon_0);

        for i = 1:num_rx
            % Ray tracing for Tx1 and Tx2
            rays_tx1 = raytrace(tx1, rx_positions{i}, pm1);
            rays_tx2 = raytrace(tx2, rx_positions{i}, pm2);
            cir_tx1 = zeros(1, max_cir_length);
            cir_tx2 = zeros(1, max_cir_length);

%% Tx1 signal
if ~isempty(rays_tx1)
    % 기존 Path Loss
    pathLoss = [rays_tx1{1, 1}.PathLoss];
    
    % 반사면 중심점 및 송신기 위치
    tx_position = [tx1.Latitude, tx1.Longitude, tx1.AntennaHeight];
    reflection_center = [37.5586876, 127.0007044, 1]; % 반사면 중심
    surface_normal = [0, 0, 1]; % 반사면 법선 벡터 (z축 기준 수직)
    
    % 입사 벡터 계산
    v_incident = reflection_center - tx_position; % 송신기 -> 반사면 중심
    k_in = v_incident / norm(v_incident); % 정규화
    
    % 반사 벡터 계산
    v_reflected = v_incident - 2 * dot(v_incident, surface_normal) * surface_normal; % 반사
    k_out = v_reflected / norm(v_reflected); % 정규화
    
    % 입사각 계산 (법선과 입사 벡터 사이의 각도)
    theta_i = acos(dot(-k_in, surface_normal)); % 입사각 (rad)
    
    % Fresnel 반사 계수 계산
    R_H = (cos(theta_i) - sqrt(complex_permittivity - sin(theta_i)^2)) / ...
          (cos(theta_i) + sqrt(complex_permittivity - sin(theta_i)^2));
    R_V = (complex_permittivity * cos(theta_i) - sqrt(complex_permittivity - sin(theta_i)^2)) / ...
          (complex_permittivity * cos(theta_i) + sqrt(complex_permittivity - sin(theta_i)^2));
    
    % 수평 및 수직 벡터 계산
    n = surface_normal / norm(surface_normal); % 법선 벡터 정규화
    H_in = cross(k_in, n) / norm(cross(k_in, n)); % 입사 방향 수평 벡터
    V_in = cross(H_in, k_in); % 입사 방향 수직 벡터
    H_rx = cross(k_out, n) / norm(cross(k_out, n)); % 반사 방향 수평 벡터
    V_rx = cross(H_rx, k_out); % 반사 방향 수직 벡터
    
    % 편광 변환 행렬 R 계산
    R = [dot(H_in, H_rx) * R_H, dot(H_in, V_rx) * R_V;
         dot(V_in, H_rx) * R_H, dot(V_in, V_rx) * R_V];
    
    % 편광 벡터 정의 (송신기 및 수신기)
    J_t = sqrt(2)/2 * [1; 1]; % 송신기 편광 벡터
    J_r = sqrt(2)/2 * [1; 1]; % 수신기 편광 벡터
    
    % 편광 손실 계산
    I_L = -20 * log10(abs(J_r' / R * J_t)); % R의 역행렬 적용
    
    % 총 Path Loss 계산
    totalPathLoss_dB = pathLoss + I_L; % 기존 경로 손실에 편광 손실 추가
    pathGains = sqrt(10.^(-totalPathLoss_dB / 10)); % dB -> 선형
    
    % 페이딩 적용
    fadingCoefficients = (randn(size(pathGains)) + 1i * randn(size(pathGains))) / sqrt(2);
    pathGains = pathGains .* fadingCoefficients;

    % CIR 업데이트
    pathDelays = [rays_tx1{1, 1}.PropagationDelay];
    for k = 1:length(pathDelays)
        [~, delayIdx] = min(abs(t - pathDelays(k)));
        if delayIdx <= length(cir_tx1)
            cir_tx1(delayIdx) = cir_tx1(delayIdx) + pathGains(k);
        end
    end
end



            %% Tx2 interference
            if ~isempty(rays_tx2)
                pathLoss = [rays_tx2{1, 1}.PathLoss];
                totalPathLoss_dB = pathLoss; % No reflection loss for Tx2
                pathGains = sqrt(10.^(-totalPathLoss_dB / 10)); % Convert to linear scale
                fadingCoefficients = (randn(size(pathGains)) + 1i * randn(size(pathGains))) / sqrt(2);
                pathGains = pathGains .* fadingCoefficients;

                pathDelays = [rays_tx2{1, 1}.PropagationDelay];
                for k = 1:length(pathDelays)
                    [~, delayIdx] = min(abs(t - pathDelays(k)));
                    if delayIdx <= length(cir_tx2)
                        cir_tx2(delayIdx) = cir_tx2(delayIdx) + pathGains(k);
                    end
                end
            end

            %% SINR Calculation
            fft_rtchan_tx1 = fft(cir_tx1);
            fft_rtchan_tx2 = fft(cir_tx2);
            signal_power = mean(abs(fft_rtchan_tx1).^2);
            interference_power = mean(abs(fft_rtchan_tx2).^2);
            sinr = 10 * log10(signal_power / (interference_power + noise_variance_watt));
            sinr_results(i, e_idx, c_idx) = sinr;
        end
    end
end

%% Plot SINR CDF
figure(1);
hold on;
for e_idx = 1:length(epsilon_r_real_vec)
    for c_idx = 1:length(conductivity_vec)
        sinr_values = squeeze(sinr_results(:, e_idx, c_idx));
        [cdf_y, cdf_x] = ecdf(sinr_values);
        plot(cdf_x, cdf_y, 'DisplayName', ...
             sprintf('\\epsilon_r = %.1f, \\sigma = %.1e', epsilon_r_real_vec(e_idx), conductivity_vec(c_idx)));
    end
end
xlabel('SINR [dB]');
ylabel('CDF');
title('SINR CDF for Different Permittivity and Conductivity');
legend show;
grid on;

figure(2);
scatter3(cellfun(@(x) x.Latitude, rx_positions), ...
         cellfun(@(x) x.Longitude, rx_positions), ...
         squeeze(mean(sinr_results, [2, 3])), 20, squeeze(mean(sinr_results, [2, 3])), 'filled');
xlabel('Latitude');
ylabel('Longitude');
zlabel('Mean SINR [dB]');
%title('Spatial Distribution of SINR');
colorbar;
grid on;
