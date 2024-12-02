clear; clc; close all;


fc = 2.4e9; 
num_rx = 1e3; 
min_radius = 10; % (m)
max_radius = 50; % (m)
txpower = 30; % Tx power: 30 dBm (1 W)
noise_psd_dbm_per_hz = -174; % dBm/Hz
path_loss_exponent = 3; 
B = 20e6; % (Hz)
noise_variance_watt = 10^((noise_psd_dbm_per_hz - 30) / 10) * B; % noise power (W)


epsilon_r_real_vec = [5, 10, 15]; % Relative permittivity for metasurfaces
conductivity_vec = [1e-4, 1e-2, 1e-1]; % Conductivity for typical RIS materials

% SINR 저장
sinr_results = zeros(num_rx, length(epsilon_r_real_vec), length(conductivity_vec));

tx1 = txsite(Name="Tx1", Latitude=37.5584876, Longitude=126.9997299, ...
    AntennaHeight=3, TransmitterPower=30, TransmitterFrequency=fc);

rx_positions = cell(1, num_rx);
max_cir_length = 0;


rx_distances = zeros(1, num_rx);
for i = 1:num_rx
    % 무작위로 반경과 각도 
    r = sqrt(rand() * (max_radius^2 - min_radius^2) + min_radius^2);
    theta = 2 * pi * rand();

    
    rx_latitude = tx1.Latitude + (r / 111000) * cos(theta);
    rx_longitude = tx1.Longitude + (r / (111000 * cosd(tx1.Latitude))) * sin(theta);
    rx_positions{i} = rxsite(Name="RX", Latitude=rx_latitude, Longitude=rx_longitude, AntennaHeight=1.7);

    
    rx_distances(i) = r;

    % propagation model
    pm1 = propagationModel("raytracing", ...
        MaxNumReflections=1, AngularSeparation="low", ...
        BuildingsMaterial="custom", ...
        BuildingsMaterialPermittivity=epsilon_r_real_vec(1), ...
        BuildingsMaterialConductivity=conductivity_vec(1), ...
        TerrainMaterial="vegetation", ...
        SurfaceMaterial="plasterboard");

    
    rays_tx1 = raytrace(tx1, rx_positions{i}, pm1);

    
    if ~isempty(rays_tx1)
        pathDelays_tx1 = [rays_tx1{1, 1}.PropagationDelay];
        max_cir_length = max(max_cir_length, ceil(max(pathDelays_tx1) * 1e9) + 1);
    end
end

t = 0:1e-9:(max_cir_length - 1) * 1e-9; 


capacity_results = zeros(num_rx, length(epsilon_r_real_vec), length(conductivity_vec));
for e_idx = 1:length(epsilon_r_real_vec)
    parfor c_idx = 1:length(conductivity_vec)
        
        epsilon_r_real = epsilon_r_real_vec(e_idx);
        conductivity = conductivity_vec(c_idx);

        pm1 = propagationModel("raytracing", ...
            MaxNumReflections=1, AngularSeparation="low", ...
            BuildingsMaterial="custom", ...
            BuildingsMaterialPermittivity=epsilon_r_real, ...
            BuildingsMaterialConductivity=conductivity, ...
            TerrainMaterial="vegetation", ...
            SurfaceMaterial="plasterboard");

        
        for i = 1:num_rx
            
            rays_tx1 = raytrace(tx1, rx_positions{i}, pm1);
            cir_tx1 = zeros(1, max_cir_length);

           
            interference_power = 0;
            num_interferers = min(10, num_rx - 1); 
            interferer_indices = randperm(num_rx, num_interferers);
            for j = interferer_indices
                if j ~= i
                    rays_interference = raytrace(tx1, rx_positions{j}, pm1);
                    cir_interference = zeros(1, max_cir_length);
                    if ~isempty(rays_interference)
                        pathDelays = [rays_interference{1, 1}.PropagationDelay];
                        pathLoss = [rays_interference{1, 1}.PathLoss];
                        
                        distance = rx_distances(j);
                        pathLoss = pathLoss + 10 * path_loss_exponent * log10(distance);
                        pathGains = sqrt(10.^(-pathLoss / 10));
                        fadingCoefficients = (randn(size(pathGains)) + 1i * randn(size(pathGains))) / sqrt(2);
                        pathGains = pathGains .* fadingCoefficients;

                        for k = 1:length(pathDelays)
                            phaseShift = exp(1j * rays_interference{1, 1}(k).PhaseShift);
                            [~, delayIdx] = min(abs(t - pathDelays(k)));
                            if delayIdx <= length(cir_interference)
                                cir_interference(delayIdx) = cir_interference(delayIdx) + pathGains(k) * phaseShift;
                            end
                        end
                    end
                    interference_power = interference_power + mean(abs(cir_interference).^2);
                end
            end

            
            if ~isempty(rays_tx1)
                pathDelays = [rays_tx1{1, 1}.PropagationDelay];
                pathLoss = [rays_tx1{1, 1}.PathLoss];
                
                distance = rx_distances(i);
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

            
            signal_power = mean(abs(cir_tx1).^2);

            % SINR 
            sinr = 10 * log10(signal_power / (interference_power + noise_variance_watt));
            sinr_results(i, e_idx, c_idx) = sinr;
        end
    end
end

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
