clear; clc; close all;

% 샘플 개수 설정
num_samples = 1e6;

% 복소 신호 생성: Rayleigh 분포를 따르기 위해 각각 독립적인 Gaussian 생성
real_part = randn(1, num_samples);
imag_part = randn(1, num_samples);
complex_signal = real_part + 1i * imag_part;

% 신호의 진폭 계산 (Rayleigh 분포를 따름)
amplitude = abs(complex_signal);

% 신호의 진폭 제곱 계산 (지수 분포를 따름)
power = amplitude .^ 2;

% 히스토그램 및 PDF 비교
figure;

% 진폭 분포 (Rayleigh 분포 확인)
subplot(2, 1, 1);
histogram(amplitude, 100, 'Normalization', 'pdf');
hold on;
x = linspace(0, max(amplitude), 1000);
rayleigh_pdf = (x ./ (1)).*exp(-(x.^2)/(2)); % Rayleigh PDF (scale parameter = 1)
plot(x, rayleigh_pdf, 'r', 'LineWidth', 1.5);
title('Rayleigh Distribution of Signal Amplitude');
xlabel('Amplitude');
ylabel('PDF');
legend('Empirical Data', 'Theoretical Rayleigh PDF');
grid on;

% 진폭 제곱 분포 (지수 분포 확인)
subplot(2, 1, 2);
histogram(power, 100, 'Normalization', 'pdf');
hold on;
x = linspace(0, max(power), 1000);
exponential_pdf = exp(-x); % Exponential PDF (mean = 1)
plot(x, exponential_pdf, 'r', 'LineWidth', 1.5);
title('Exponential Distribution of Squared Amplitude');
xlabel('Power (Amplitude^2)');
ylabel('PDF');
legend('Empirical Data', 'Theoretical Exponential PDF');
grid on;
