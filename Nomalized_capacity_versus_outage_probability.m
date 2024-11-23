clear; clc; close all;

fc = 2.4e9; % 송신 주파수: 2.4 GHz
num_rx = 1e3; % 각 반지름 내 배치할 RX 개수
min_radius = 10; % 최소 반지름 (m)
max_radius = 50; % 최대 반지름 (m)
txpower = 30; % Tx 파워: 30 dBm (1 W)
noise_psd_dbm_per_hz = -174; % dBm/Hz
path_loss_exponent = 3; % path loss exponent
B = 20e6; % Bandwidth (Hz)
noise_variance_watt = 10^((noise_psd_dbm_per_hz - 30) / 10) * B; % noise power (W)

% RIS reflection coefficient 설정 
alpha_vec=[0.1 0.3 0.5 0.7 0.9];
sinr=[];

%osm_file_path = 'C:\Users\kimyj\OneDrive - 동국대학교\바탕 화면\raytracing_jhcodedongguk.osm'; 
%sv = siteviewer("Buildings", "dongguk.osm");

% 송신기 1 (Tx1) 설정 (left)
tx1 = txsite(Name="Tx1", ...
    Latitude=37.5584876, ...
    Longitude=126.9997299, ... %127.0001671 %left: decreasing longitude
    AntennaHeight=3, ...
    TransmitterPower=txpower, ...
    TransmitterFrequency=fc);
%show(tx1, "ShowAntennaHeight", true);

% 송신기 2 (Tx2 - 간섭원) 설정
tx2 = txsite(Name="Tx2", ...
    Latitude=37.5582876, ...
    Longitude=127.0009671, ...
    AntennaHeight=3, ...
    TransmitterPower=txpower, ...
    TransmitterFrequency=fc);
%show(tx2, "ShowAntennaHeight", true);

% pm1 설정 (RIS 반사 경로)
pm1 = propagationModel("raytracing", ...
    MaxNumReflections=2, ... 
    AngularSeparation="medium", ...
    BuildingsMaterial="perfect-reflector");

% pm2 설정 (간섭 경로)
pm2 = propagationModel("raytracing", ...
    MaxNumReflections=2, ...
    AngularSeparation="medium", ...
    BuildingsMaterial="perfect-reflector");

% RX 위치 설정 및 분석
rx_positions = cell(1, num_rx); % RX 위치 객체 저장
sinr_values = zeros(1, num_rx); % SINR 값 저장
max_cir_length = 0; % CIR의 최대 길이 초기화

for g=1:length(alpha_vec)

    alpha=alpha_vec(g);

% RX를 최소 반경과 최대 반경 사이에 무작위로 균등하게 배치
for i = 1:num_rx
    % 무작위로 반지름과 각도 생성 (도넛 모양 배치)
    r = sqrt(rand() * (max_radius^2 - min_radius^2) + min_radius^2); % 무작위 반경 생성
    theta = 2 * pi * rand(); % 무작위 각도 생성
    
    % RX의 위도 및 경도 계산
    rx_latitude = tx1.Latitude + (r / 111000) * cos(theta);
    rx_longitude = tx1.Longitude + (r / (111000 * cosd(tx1.Latitude))) * sin(theta);
    rx_positions{i} = rxsite(Name="RX", Latitude=rx_latitude, Longitude=rx_longitude, AntennaHeight=1.7);
    %show(rx, "ShowAntennaHeight", true);

    % CIR 계산 준비
    rays_tx1 = raytrace(tx1, rx_positions{i}, pm1); % pm1을 통해 Tx1의 신호 경로
    rays_tx2 = raytrace(tx2, rx_positions{i}, pm2); % pm2를 통해 Tx2의 신호 경로 (간섭원)

    % 최대 CIR 길이 갱신
    if ~isempty(rays_tx1)
        pathDelays = [rays_tx1{1, 1}.PropagationDelay];
        fs = 1e9; %channel resolution issue 
        t = 0:1/fs:max(pathDelays) + 1/fs;
        cir_length = length(t);
        if cir_length > max_cir_length
            max_cir_length = cir_length;
        end
    end
end

% CIR 계산 및 SINR 분석
for i = 1:num_rx
    % Tx1 신호 분석
    rays_tx1 = raytrace(tx1, rx_positions{i}, pm1);
    rays_tx1 = rays_tx1{1, 1};
    cir_tx1 = zeros(1, max_cir_length);

    if ~isempty(rays_tx1)
        pathDelays = [rays_tx1.PropagationDelay];
        pathLoss = [rays_tx1.PathLoss];
        
        % 거리 기반 path loss exponent 적용
        distance = norm([tx1.Latitude - rx_positions{i}.Latitude, tx1.Longitude - rx_positions{i}.Longitude]) * 111000;
        pathLoss = pathLoss + 10 * path_loss_exponent * log10(distance);

        pathGains = sqrt(10.^(-pathLoss / 10));
        fadingCoefficients = (randn(size(pathGains)) + 1i * randn(size(pathGains))) / sqrt(2);
        pathGains = pathGains .* fadingCoefficients;

        % CIR 계산 (RIS 반사계수 적용)
        fs = 1e9;
        t = 0:1/fs:(max_cir_length - 1) / fs;
        for k = 1:length(pathDelays)
            [~, delayIdx] = min(abs(t - pathDelays(k)));
            if delayIdx <= length(cir_tx1)
                % 반사계수 alpha 적용
                cir_tx1(delayIdx) = cir_tx1(delayIdx) + alpha * pathGains(k);
            end
        end
    end

    % Tx2 신호 (간섭 신호) 분석
    rays_tx2 = raytrace(tx2, rx_positions{i}, pm2);
    rays_tx2 = rays_tx2{1, 1};
    cir_tx2 = zeros(1, max_cir_length);

    if ~isempty(rays_tx2)
        pathDelays = [rays_tx2.PropagationDelay];
        pathLoss = [rays_tx2.PathLoss];
        
        % 거리 기반 path loss exponent 적용
        distance = norm([tx2.Latitude - rx_positions{i}.Latitude, tx2.Longitude - rx_positions{i}.Longitude]) * 111000;
        pathLoss = pathLoss + 10 * path_loss_exponent * log10(distance);

        pathGains = sqrt(10.^(-pathLoss / 10));
        fadingCoefficients = (randn(size(pathGains)) + 1i * randn(size(pathGains))) / sqrt(2);
        pathGains = pathGains .* fadingCoefficients;

        % CIR 계산
        for k = 1:length(pathDelays)
            [~, delayIdx] = min(abs(t - pathDelays(k)));
            if delayIdx <= length(cir_tx2)
                cir_tx2(delayIdx) = cir_tx2(delayIdx) + pathGains(k);
            end
        end
    end

    % Tx1으로부터의 신호 파워 (신호 파워)
    fft_rtchan_tx1 = fft(cir_tx1);
    Y_signal = ones(1, length(fft_rtchan_tx1)) .* fft_rtchan_tx1;
    signal_power = mean(abs(Y_signal).^2);

    % Tx2로부터의 신호 파워 (간섭 파워)
    fft_rtchan_tx2 = fft(cir_tx2);
    Y_interference = ones(1, length(fft_rtchan_tx2)) .* fft_rtchan_tx2;
    interference_power = mean(abs(Y_interference).^2);

    % SINR 계산
    sinr_values(i) = 10 * log10(signal_power / (interference_power + noise_variance_watt));
end

sinr=[sinr sinr_values.'];

end

%% Outage Probability와 Normalized Capacity (C/B) 그래프 그리기
figure;
hold on;
set(gca, 'XScale', 'log');
xlabel('Outage Probability');
ylabel('Normalized Capacity (C/B)');
title('Nomrmalized Capacity versus Outage Probability');
box on;
grid on;

% 각 반사 계수에 대해 그래프 그리기
for g = 1:length(alpha_vec)
    sinr_values = sinr(:, g);

    % 실제 Outage Probability 계산: 임계값 이하인 경우의 비율을 계산하여 여러 임계값에 대해 반복
    outage_prob = [];
    capacity = [];

    % 임계값 범위 설정 (SINR 최소값부터 최대값까지 균일하게 분포)
    threshold_db = linspace(min(sinr_values), max(sinr_values), 100);

    for k = 1:length(threshold_db)
        % 현재 임계값에서의 Outage Probability 계산
        current_outage_prob = sum(sinr_values < threshold_db(k)) / num_rx;
        outage_prob = [outage_prob, current_outage_prob];

        % 현재 임계값에서의 용량 계산 (정규화된 용량)
        % 이때 임계값 threshold_db(k)를 사용하여 용량 계산
        current_capacity = log2(1 + 10^(threshold_db(k) / 10));
        capacity = [capacity, current_capacity];
    end

    % 그래프 그리기
    plot(outage_prob, capacity, 'LineWidth', 1.5);
end

legend(arrayfun(@(x) sprintf('Reflection coefficient=%.1f', x), alpha_vec, 'UniformOutput', false), 'Location', 'best');
