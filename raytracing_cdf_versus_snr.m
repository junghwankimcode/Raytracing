clear; clc; close all;

% 기본 설정
fc = 2.4e9; % 송신 주파수: 2.4 GHz
num_rx = 1e3; % 각 반지름 내 배치할 RX 개수
min_radius = 10; % 최소 반지름 (m)
max_radius = 50; % 최대 반지름 (m)
txpower = 30; % Tx 파워: 30 dBm (1 W)
noise_psd_dbm_per_hz = -174; % dBm/Hz
path_loss_exponent = 3; % path loss exponent
B = 60e6; % Bandwidth (Hz)
noise_variance_watt = 10^((noise_psd_dbm_per_hz - 30) / 10) * B; % noise power (Watt)
% noise_variance_watt = noise_psd_dbm_per_hz + 10*log10(60e6);


% Site Viewer 설정
%osm_file_path = 'C:\Users\kimyj\OneDrive - 동국대학교\바탕 화면\raytracing_jhcodedongguk.osm'; % OSM 파일의 경로
%sv = siteviewer("Buildings", "dongguk.osm");

% 송신기 설정
tx = txsite(Name="Standard transmitter", ...
    Latitude=37.5581876, ...
    Longitude=127.0001671, ...
    AntennaHeight=3, ...
    TransmitterPower=txpower, ...
    TransmitterFrequency=fc);
%show(tx, "ShowAntennaHeight", true);
pm = propagationModel("raytracing", ...
    MaxNumReflections=4, ...
    MaxNumDiffractions=2,...
    AngularSeparation="medium", ...
    BuildingsMaterial="concrete", ...
    TerrainMaterial="vegetation");


% RX 위치 설정 및 분석
rx_positions = cell(1, num_rx); % RX 위치 객체 저장
actualSNR_values = zeros(1, num_rx); % 실제 SNR 값 저장
max_cir_length = 0; % CIR의 최대 길이 초기화

% RX를 최소 반경과 최대 반경 사이에 무작위로 균등하게 배치
for i = 1:num_rx
    % 무작위로 반지름과 각도 생성 (도넛 모양 배치)
    r = sqrt(rand() * (max_radius^2 - min_radius^2) + min_radius^2); % 무작위 반경 생성
    theta = 2 * pi * rand(); % 무작위 각도 생성
    
    % RX의 위도 및 경도 계산
    rx_latitude = tx.Latitude + (r / 111000) * cos(theta);
    rx_longitude = tx.Longitude + (r / (111000 * cosd(tx.Latitude))) * sin(theta);
    rx_positions{i} = rxsite(Name="RX", Latitude=rx_latitude, Longitude=rx_longitude, AntennaHeight=1.7);
    %show(rx_positions{i}, "ShowAntennaHeight", true);
    % CIR 계산 준비
    rays = raytrace(tx, rx_positions{i}, pm);
    rays = rays{1, 1};
    if ~isempty(rays)
        pathDelays = [rays.PropagationDelay];
        fs = 1e9;
        t = 0:1/fs:max(pathDelays) + 1/fs;
        cir_length = length(t);

        % 최대 CIR 길이 갱신
        if cir_length > max_cir_length
            max_cir_length = cir_length;
        end
    end
end

% CIR 계산 및 분석
for i = 1:num_rx
    rays = raytrace(tx, rx_positions{i}, pm);
    rays = rays{1, 1};
    if ~isempty(rays)
        pathDelays = [rays.PropagationDelay];
        pathLoss = [rays.PathLoss];

        % path loss exponent 적용 (dB 스케일에서 거리 기반으로 적용)
        distance = norm([tx.Latitude - rx_positions{i}.Latitude, tx.Longitude - rx_positions{i}.Longitude]) * 111000; % 거리 계산 (미터 단위)
        pathLoss = pathLoss + 10 * path_loss_exponent * log10(distance);

        pathGains = sqrt(10.^(-pathLoss / 10));

        % 페이딩 계수 추가 (Rayleigh 페이딩)
        fadingCoefficients = (randn(size(pathGains)) + 1i * randn(size(pathGains))) / sqrt(2);
        pathGains = pathGains .* fadingCoefficients;

        % CIR 계산
        fs = 1e9;
        t = 0:1/fs:(max_cir_length - 1) / fs;
        cir = zeros(1, max_cir_length);

        for k = 1:length(pathDelays)
            [~, delayIdx] = min(abs(t - pathDelays(k)));
            if delayIdx <= length(cir)
                cir(delayIdx) = cir(delayIdx) + pathGains(k);
            end
        end

        % FFT 기반 주파수 응답 계산
        fft_rtchan = fft(cir);
        Y = ones(1, length(fft_rtchan)) .* fft_rtchan; % 주파수 응답 신호
        
        noise = sqrt(noise_variance_watt / 2) * (randn(size(Y)) + 1i * randn(size(Y)));

        %delta_f = fs / length(Y);

        % 실제 SNR 계산
        signal_power = mean(abs(Y).^2);
        actualSNR = 10 * log10(signal_power / noise_variance_watt);
        actualSNR_values(i) = actualSNR;
    else
        actualSNR_values(i) = NaN;
    end
end

% NaN 제거
actualSNR_values = actualSNR_values(~isnan(actualSNR_values));

%% SNR 분포 히스토그램
figure;
histogram(actualSNR_values, 'Normalization', 'probability', 'BinWidth', 0.65);
xlabel('SNR [dB]');
ylabel('Empirical probability density');
%title('PDF of SNR');
grid on;

% SNR CDF 계산 및 그래프
[sortedSNR, ~] = sort(actualSNR_values);
cdf_values = (1:length(sortedSNR)) / length(sortedSNR);
figure;
plot(sortedSNR, cdf_values, 'LineWidth', 1.5);
xlabel('SNR [dB]');
ylabel('Cumulative probability');
%title('CDF of SNR');
grid on;

% 25% 및 75% 백분위수 계산
snr_25 = prctile(actualSNR_values, 25);
snr_75 = prctile(actualSNR_values, 75);

% 25% 및 75% 위치에 마커 추가
hold on;
plot(snr_25, 0.25, 'r*', 'MarkerSize', 8);
plot(snr_75, 0.75, 'bd', 'MarkerSize', 8);
legend('25th Percentile','75th Percentile');

% 25%와 75% SNR 값 출력
disp(['25th Percentile SNR: ', num2str(snr_25), ' dB']);
disp(['75th Percentile SNR: ', num2str(snr_75), ' dB']);
