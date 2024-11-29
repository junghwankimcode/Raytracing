clear; clc; close all;

% 파라미터 설정
fc = 2.4e9; % 송신 주파수: 2.4 GHz
num_rx = 1e3; % RX 개수
min_radius = 10; % 최소 반지름 (m)
max_radius = 50; % 최대 반지름 (m)
txpower = 30; % Tx 파워: 30 dBm (1 W)
noise_psd_dbm_per_hz = -174; % dBm/Hz
path_loss_exponent = 3; % 경로 손실 지수
B = 20e6; % 대역폭 (Hz)
noise_variance_watt = 10^((noise_psd_dbm_per_hz - 30) / 10) * B; % noise power (W)

% 유전율 및 전도도 설정
epsilon_r_real_vec = [5, 10, 15]; % Relative permittivity for metasurfaces
conductivity_vec = [1e-4, 1e-2, 1e-1]; % Conductivity for typical RIS materials


% SINR 저장
sinr_results = zeros(num_rx, length(epsilon_r_real_vec), length(conductivity_vec));

% 송신기 설정
tx1 = txsite(Name="Tx1", Latitude=37.5584876, Longitude=126.9997299, ...
    AntennaHeight=3, TransmitterPower=30, TransmitterFrequency=fc);

tx2 = txsite(Name="Tx2", Latitude=37.5582876, Longitude=127.0009671, ...
    AntennaHeight=3, TransmitterPower=30, TransmitterFrequency=fc);

% RX 위치 설정
rx_positions = cell(1, num_rx);
max_cir_length = 0;

% RX 위치 초기화
for i = 1:num_rx
    % 무작위로 반경과 각도 생성
    r = sqrt(rand() * (max_radius^2 - min_radius^2) + min_radius^2);
    theta = 2 * pi * rand();

    % RX 위치 계산
    rx_latitude = tx1.Latitude + (r / 111000) * cos(theta);
    rx_longitude = tx1.Longitude + (r / (111000 * cosd(tx1.Latitude))) * sin(theta);
    rx_positions{i} = rxsite(Name="RX", Latitude=rx_latitude, Longitude=rx_longitude, AntennaHeight=1.7);

    % CIR 길이 계산 (Tx1과 Tx2 모두 고려)
    rays_tx1 = raytrace(tx1, rx_positions{i}, propagationModel("raytracing", MaxNumReflections=2));
    rays_tx2 = raytrace(tx2, rx_positions{i}, propagationModel("raytracing", MaxNumReflections=2));

    % CIR 길이 갱신
    if ~isempty(rays_tx1)
        pathDelays_tx1 = [rays_tx1{1, 1}.PropagationDelay];
        max_cir_length = max(max_cir_length, ceil(max(pathDelays_tx1) * 1e9) + 1);
    end
    if ~isempty(rays_tx2)
        pathDelays_tx2 = [rays_tx2{1, 1}.PropagationDelay];
        max_cir_length = max(max_cir_length, ceil(max(pathDelays_tx2) * 1e9) + 1);
    end
end

t = 0:1e-9:(max_cir_length - 1) * 1e-9; % 시간 벡터 생성

% SINR 및 용량 계산
capacity_results = zeros(num_rx, length(epsilon_r_real_vec), length(conductivity_vec));
for e_idx = 1:length(epsilon_r_real_vec)
    for c_idx = 1:length(conductivity_vec)
        % 사용자 정의 재질 설정
        epsilon_r_real = epsilon_r_real_vec(e_idx);
        conductivity = conductivity_vec(c_idx);

        pm1 = propagationModel("raytracing", ...
            MaxNumReflections=2, AngularSeparation="medium", ...
            BuildingsMaterial="custom", ...
            BuildingsMaterialPermittivity=epsilon_r_real, ...
            BuildingsMaterialConductivity=conductivity,...
            TerrainMaterial="vegetation",...
            SurfaceMaterial="plasterboard");

        pm2 = propagationModel("raytracing", ...
            MaxNumReflections=2, AngularSeparation="medium", ...
            BuildingsMaterial="concrete",...
            TerrainMaterial="vegetation",...
            SurfaceMaterial="plasterboard");

        % 각 RX 위치에서 SINR 및 용량 계산
        for i = 1:num_rx
            % Tx1 (RIS) 레이 트레이싱
            rays_tx1 = raytrace(tx1, rx_positions{i}, pm1);
            cir_tx1 = zeros(1, max_cir_length);

            % Tx2 (간섭원) 레이 트레이싱
            rays_tx2 = raytrace(tx2, rx_positions{i}, pm2);
            cir_tx2 = zeros(1, max_cir_length);

            %% CIR 계산
            % Tx1 (RIS)
            if ~isempty(rays_tx1)
                pathDelays = [rays_tx1{1, 1}.PropagationDelay];
                pathLoss = [rays_tx1{1, 1}.PathLoss];
                % 거리 기반 path loss exponent 적용
                distance = norm([tx1.Latitude - rx_positions{i}.Latitude, tx1.Longitude - rx_positions{i}.Longitude]) * 111000;
                pathLoss = pathLoss + 10 * path_loss_exponent * log10(distance);
        
                pathGains = sqrt(10.^(-pathLoss / 10));
                fadingCoefficients = (randn(size(pathGains)) + 1i * randn(size(pathGains))) / sqrt(2);
                pathGains = pathGains .* fadingCoefficients;

                for k = 1:length(pathDelays)
                    phaseShift = exp(1j * rays_tx1{1, 1}(k).PhaseShift);
                    [~, delayIdx] = min(abs(t - pathDelays(k)));
                    if delayIdx <= length(cir_tx1)
                        cir_tx1(delayIdx) = cir_tx1(delayIdx) + pathGains(k) * phaseShift;
                    end
                end
            end

            % Tx2 (간섭원)
            if ~isempty(rays_tx2)
                pathDelays = [rays_tx2{1, 1}.PropagationDelay];
                pathLoss = [rays_tx2{1, 1}.PathLoss];
                % 거리 기반 path loss exponent 적용
                distance = norm([tx2.Latitude - rx_positions{i}.Latitude, tx2.Longitude - rx_positions{i}.Longitude]) * 111000;
                pathLoss = pathLoss + 10 * path_loss_exponent * log10(distance);
        
                pathGains = sqrt(10.^(-pathLoss / 10));
                fadingCoefficients = (randn(size(pathGains)) + 1i * randn(size(pathGains))) / sqrt(2);
                pathGains = pathGains .* fadingCoefficients;

                for k = 1:length(pathDelays)
                    phaseShift = exp(1j * rays_tx2{1, 1}(k).PhaseShift);
                    [~, delayIdx] = min(abs(t - pathDelays(k)));
                    if delayIdx <= length(cir_tx2)
                        cir_tx2(delayIdx) = cir_tx2(delayIdx) + pathGains(k) * phaseShift;
                    end
                end
            end

            %% 신호 및 간섭 파워 계산
            signal_power = mean(abs(cir_tx1).^2);
            interference_power = mean(abs(cir_tx2).^2);

            % SINR 계산
            sinr = 10 * log10(signal_power / (interference_power + noise_variance_watt));
            sinr_results(i, e_idx, c_idx) = sinr;

            % 용량 계산
            capacity_results(i, e_idx, c_idx) = B * log2(1 + 10^(sinr / 10));
        end
    end
end


% %% 유전율 기준 SINR 그래프
% figure;
% hold on;
% mean_sinr_per_epsilon = squeeze(mean(sinr_results, [1, 3])); % 전도도 축 평균
% plot(epsilon_r_real_vec, mean_sinr_per_epsilon, '-o', 'LineWidth', 1.5);
% xlabel('Relative Permittivity (\epsilon_r)');
% ylabel('Mean SINR [dB]');
% title('Mean SINR vs Relative Permittivity (\epsilon_r)');
% grid on;
% 
% %% 전도도 기준 SINR 그래프
% figure;
% hold on;
% mean_sinr_per_sigma = squeeze(mean(sinr_results, [1, 2])); % 유전율 축 평균
% plot(conductivity_vec, mean_sinr_per_sigma, '-s', 'LineWidth', 1.5);
% xlabel('Conductivity (\sigma) [S/m]');
% ylabel('Mean SINR [dB]');
% title('Mean SINR vs Conductivity (\sigma)');
% grid on;



%%
% [X, Y] = meshgrid(epsilon_r_real_vec, conductivity_vec);
% Z = squeeze(mean(capacity_results, 1)); % RX 평균 Capacity
% 
% figure;
% surf(X, Y, Z, 'EdgeColor', 'none');
% xlabel('Relative Permittivity (\epsilon_r)');
% ylabel('Conductivity (\sigma) [S/m]');
% zlabel('Capacity (bps)');
% % title('3D Plot: Capacity vs \epsilon_r and \sigma');
% colorbar;
% 
% 
% %%
% [X, Y] = meshgrid(epsilon_r_real_vec, conductivity_vec);
% Z = squeeze(mean(sinr_results, 1)); % RX 평균 Capacity
% 
% figure;
% surf(X, Y, Z, 'EdgeColor', 'none');
% xlabel('Relative Permittivity (\epsilon_r)');
% ylabel('Conductivity (\sigma) [S/m]');
% zlabel('SINR [dB]');
% % title('3D Plot: Capacity vs \epsilon_r and \sigma');
% colorbar;

figure;
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

% RX 위치별 SINR 히트맵
figure;
scatter3(cellfun(@(x) x.Latitude, rx_positions), ...
         cellfun(@(x) x.Longitude, rx_positions), ...
         squeeze(mean(sinr_results, [2, 3])), 20, squeeze(mean(sinr_results, [2, 3])), 'filled');
xlabel('Latitude');
ylabel('Longitude');
zlabel('Mean SINR [dB]');
title('Spatial Distribution of SINR');
colorbar;
grid on;

