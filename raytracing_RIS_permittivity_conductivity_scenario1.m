clear; clc; close all;

fc = 2.4e9; 
num_rx = 1e3; 
min_radius = 10; % (m)
max_radius = 50; % (m)
txpower = 30; % 30 dBm (1 W)
noise_psd_dbm_per_hz = -174; % dBm/Hz
path_loss_exponent = 3; 
B = 20e6; % (Hz)
noise_variance_watt = 10^((noise_psd_dbm_per_hz - 30) / 10) * B; % noise power (W)

epsilon_r_real_vec = [5, 10, 15]; % Relative permittivity for metasurfaces
conductivity_vec = [1e-4, 1e-2, 1e-1]; % Conductivity for typical RIS materials


sinr_results = zeros(num_rx, length(epsilon_r_real_vec), length(conductivity_vec));

tx1 = txsite(Name="Tx1", Latitude=37.5584876, Longitude=126.9997299, ...
    AntennaHeight=3, TransmitterPower=30, TransmitterFrequency=fc);

tx2 = txsite(Name="Tx2", Latitude=37.5582876, Longitude=127.0009671, ...
    AntennaHeight=3, TransmitterPower=30, TransmitterFrequency=fc);

rx_positions = cell(1, num_rx);
max_cir_length = 0;


for i = 1:num_rx
    
    r = sqrt(rand() * (max_radius^2 - min_radius^2) + min_radius^2);
    theta = 2 * pi * rand();

   
    rx_latitude = tx1.Latitude + (r / 111000) * cos(theta);
    rx_longitude = tx1.Longitude + (r / (111000 * cosd(tx1.Latitude))) * sin(theta);
    rx_positions{i} = rxsite(Name="RX", Latitude=rx_latitude, Longitude=rx_longitude, AntennaHeight=1.7);

   
    rays_tx1 = raytrace(tx1, rx_positions{i}, propagationModel("raytracing", MaxNumReflections=2));
    rays_tx2 = raytrace(tx2, rx_positions{i}, propagationModel("raytracing", MaxNumReflections=2));

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


capacity_results = zeros(num_rx, length(epsilon_r_real_vec), length(conductivity_vec));
for e_idx = 1:length(epsilon_r_real_vec)
    for c_idx = 1:length(conductivity_vec)
     
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

  
        for i = 1:num_rx
            
            rays_tx1 = raytrace(tx1, rx_positions{i}, pm1);
            cir_tx1 = zeros(1, max_cir_length);

            
            rays_tx2 = raytrace(tx2, rx_positions{i}, pm2);
            cir_tx2 = zeros(1, max_cir_length);

            %% CIR 계산
            % Tx1 (RIS)
            if ~isempty(rays_tx1)
                pathDelays = [rays_tx1{1, 1}.PropagationDelay];
                pathLoss = [rays_tx1{1, 1}.PathLoss];
             
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

          
            signal_power = mean(abs(cir_tx1).^2);
            interference_power = mean(abs(cir_tx2).^2);

     
            sinr = 10 * log10(signal_power / (interference_power + noise_variance_watt));
            sinr_results(i, e_idx, c_idx) = sinr;

            capacity_results(i, e_idx, c_idx) = B * log2(1 + 10^(sinr / 10));
        end
    end
end


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

