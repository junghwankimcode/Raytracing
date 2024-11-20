clear; clc; close all;


fc = 2.4e9; % 송신 주파수: 2.4 GHz
initial_distance = 10; % 첫 번째 RX와 TX 사이 거리 (미터)
distance_step = 20; % RX 간의 간격 (미터)
num_positions = 1; % RX 위치 개수
txpower = 30; % 송신기 파워: 1 W (30 dBm)
noise_psd_dbm_per_hz = -174; % dBm/Hz
path_loss_exponent = 2; % path loss exponent
B = 60e6; % Bandwidth (Hz)
fixed_noise_variance = 10^((noise_psd_dbm_per_hz - 30) / 10) * B; % noise power (Watt)

% Site Viewer 설정
osm_file_path = 'C:\Users\kimyj\OneDrive - 동국대학교\바탕 화면\raytracing\dongguk.osm'; % OSM 파일의 경로
%sv = siteviewer("Buildings", "hongkong.osm");
sv = siteviewer("Buildings", "dongguk.osm");
% 송신기 설정
% tx = txsite(Name="Standard transmitter", ...
%     Latitude=22.2789, ...
%     Longitude=114.1625, ...
%     AntennaHeight=5, ...
%     TransmitterPower=txpower, ...
%     TransmitterFrequency=fc);
tx = txsite(Name="Standard transmitter", ...
    Latitude= 37.5582876, ...
    Longitude=  127.0001671, ...
    AntennaHeight=3, ...
    TransmitterPower=txpower, ...
    TransmitterFrequency=fc);

%show(tx, "ShowAntennaHeight", true);

% 각 RX 위치 설정 및 분석
rx_positions = cell(num_positions, 1); % RX 위치 객체 저장
estiSNR_emil_values = zeros(1, num_positions); % SNR 추정 결과 저장
actualSNR_values = zeros(1, num_positions); % 실제 SNR 값 저장
Y_noisy_vec = cell(1, num_positions); % y_noisy 저장
noise_variance_vec = cell(1, num_positions);
noise_vec = cell(1, num_positions);
distance_m_vec = zeros(1, num_positions); % 실제 거리 벡터

% 고정된 노이즈 분산 값 설정 (Watt 단위)
%fixed_noise_variance = 1e-9; 

for idx = 1:num_positions
    % RX 위치 설정 (TX에서 일정 거리만큼 떨어진 직선상에 배치)
    actual_distance_m = initial_distance + (idx-1) * distance_step; % 실제 직선 거리 (미터)
    new_latitude = tx.Latitude + (actual_distance_m / 111000); % 위도 방향으로 거리 변화
    rx = rxsite(Name="Standard receiver", Latitude=new_latitude, Longitude=tx.Longitude, AntennaHeight=1.7); 
    rx_positions{idx} = rx;
    %show(rx, "ShowAntennaHeight", true);

    % 경로 모델 및 경로 추적
    pm = propagationModel("raytracing", MaxNumReflections=4, AngularSeparation="medium", BuildingsMaterial="concrete", TerrainMaterial="vegetation");
    rays = raytrace(tx, rx, pm);
    coverage(tx, pm, ...
    SignalStrengths=-40:5:30, ...  % 신호 강도 범위 (dBm)
    MaxRange=50, ...               % 최대 거리 (예: 수신기 최대 거리보다 조금 여유있게)
    Colormap='jet',...
    Resolution=5, ...              % 해상도 설정
    Transparency=0.7);             % 투명도 설정
    hide(tx);
    
    % 경로 추적 결과 시각화
    rays = rays{1, 1}; % 1번째 Tx에서 1번째 Rx까지의 경로
    %plot(rays);

    rtchan = comm.RayTracingChannel(rays, tx, rx);
    rtchan.SampleRate = 60e6;
    showProfile(rtchan)
    % 경로별 지연 및 감쇠 정보
    pathDelays = [rays.PropagationDelay]; % path delay (초)
    pathLoss = [rays.PathLoss]; % path loss (dB)

    % dB에서 선형 감쇠로 변환
    
    distance=norm([tx.Latitude-rx.Latitude,tx.Longitude-rx.Longitude])*111000;
    pathLoss=pathLoss+10*path_loss_exponent*log10(distance);

    pathGains = sqrt(10.^(-pathLoss / 10));

    fadingCoefficients = (randn(size(pathGains))+1i*randn(size(pathGains)))/sqrt(2);
    pathGains=pathGains.*fadingCoefficients;

    % 실제 거리 저장
    distance_m_vec(idx) = actual_distance_m;

    % 경로 지연 시간 계산
    c = 3e8; % 빛의 속도 (m/s)
    tau = actual_distance_m / c; % (초)

    fs = B*3; % 샘플링 주파수
    t = 0:1/fs:max(pathDelays) + 1/fs; 

    cir = zeros(size(t));
    for k = 1:length(pathDelays)
        % 각 지연 시간에 가장 가까운 인덱스를 찾기
        [~, delayIdx] = min(abs(t - pathDelays(k)));
        if delayIdx <= length(cir)
            cir(delayIdx) = cir(delayIdx) + pathGains(k); % path gain을 해당 인덱스에 추가
        end
    end

    % FFT 기반 주파수 응답
    fft_rtchan = fft(cir); 
    N = length(cir);
    f = (0:N-1) * (fs / N); 
    X_frequency = ones(1, size(fft_rtchan, 2));
    Y = X_frequency .* fft_rtchan; % original signal

    % 실제 SNR 계산
    signal_power = mean(abs(Y).^2); % 신호 전력 (Watt 단위)
    actualSNR = 10 * log10(signal_power / fixed_noise_variance); % 실제 SNR (dB 단위)
    actualSNR_values(idx) = actualSNR;

    % SNR 추정 계산
    noise = sqrt(fixed_noise_variance / 2) * (randn(size(Y)) + 1i * randn(size(Y))); % 복소 가우시안 잡음 생성
    Y_noisy = Y + noise; % 노이즈가 추가된 신호

%     trainDoF = numel(rays);
%     [~, S2, ~, ~] = makeHankel(Y_noisy);
%     S2 = diag(S2);
%     a = sum(S2(1:trainDoF).^(2)) - (trainDoF / (length(S2) - trainDoF)) * sum(S2(trainDoF+1:end).^(2));
%     b = (length(S2) / (length(S2) - trainDoF)) * sum(S2(trainDoF+1:end).^(2));
%     estiSNR_emil = 10 * log10(a / b);
%     estiSNR_emil_values(idx) = estiSNR_emil;

    Y_noisy_vec{idx} = Y_noisy;
    noise_vec{idx} = noise;
    noise_variance_vec{idx} = fixed_noise_variance;

    % 그래프 표시
    figure(idx);
    subplot(2,1,1)
    stem(t, abs(cir));
    xlabel('Delay (s)');
    ylabel('Amplitude');
    title(['Time-Domain Channel Impulse Response (CIR) - RX Position: ', num2str(idx)]);

    subplot(2,1,2)
    fft_rtchan = fft(cir); 
    N = length(cir);
    f = (0:N-1) * (fs / N); 
    plot(f, abs(fft_rtchan));
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    title(['Frequency-Domain Channel Response (CFR) - RX Position: ', num2str(idx)]);
end
