%%
% Script: noiseAnalysis_Chambers.m
%
% This script evaluates how Gaussian noise in the wheel-speed signal
% affects computed distance and regenerative energy. It uses a baseline
% (masked true velocity) to compare against noisy signals with various
% sigma levels. The script computes mean and standard deviation of error
% for distance and energy, determines the maximum allowable noise level,
% and visualizes the results with shaded error regions. The functional
% result that this script gives is the overall allowable sigma given the
% client specified tolerance that satisfies both conditions given modelled
% noise that may be experienced. 
%
% The script calls the following functions:
%
% [D, counts] = betterDifferentiator(t, v)
%   Computes numerical derivative (acceleration) from velocity data.
%   Returns derivative vector and counts of finite difference method used.
%
% [I, counts] = betterIntegrator(t, x)
%   Computes numerical integral of a signal over time.
%   Returns integral value and counts of integration method used.
% 
% This script additionally utilizes the WheelNoiseApp: 
%   Seed: 42
%   Gaussian: 0 (manually added later)
%   Quantization: .05
%   Dropout frac: .02
%   Zooms: 0
%   This app serves to generate a masked dataset given velocity and time
%   data from the UDDS test, and gives us the basis for our noise analysis.
%
% Outputs:
% - Figures showing distance and energy error vs Gaussian noise
% - Console output of baseline and max allowable sigma
%
% Author:       Isaac Chambers
% Section:      ME 2016 - C
% Assignment:   MP4
% Date:         11/23/2025

clc;
format longG;

%% Load wheel-speed data
% Ensure WheelNoiseApp has been run and table S exists in workspace
% Contains time, true velocity, and wheel velocity with dropout/quantization
if ~exist('S','var')
    error('Table S not found. Run WheelNoiseApp and send S to workspace.');
end

% Grab WheelNoiseApp generated data and convert wheel speed to m/s
t_noisy = S.time;
v_noisy = (1609.344/3600) * S.v_wheel;  % wheel signal with noise/dropout applied
v_true  = (1609.344/3600) * S.v_true;   % true velocity signal

%% Baseline calculation using masked true velocity
% Set vehicle and regen parameters
m     = 1500; 
eta   = 0.75;
Pclip = 80000;

% Compute acceleration from velocity and regen power
[accel_base, ~] = betterDifferentiator(t_noisy, v_true);
P_base          = eta .* min(Pclip, max(0, -m.*accel_base.*v_true));

% Integrate power and velocity to get energy and distance
[E_base, ~] = betterIntegrator(t_noisy, P_base);
[D_base, ~] = betterIntegrator(t_noisy, v_true);

% Convert energy to kWh
E_base_kWh   = E_base / 3.6e6;

% Print baseline results
fprintf('Baseline (masked true) results:\n');
fprintf('Distance traveled: %.3f km\n', D_base/1000);
fprintf('Regen energy recovered: %.4f kWh\n\n', E_base_kWh);

%% Noise sweep over Gaussian sigma
% Vary sigma, compute distance and energy errors, track mean/std
sigma_vals = linspace(0, .25, 500); 
dist_error_mean   = zeros(size(sigma_vals));
dist_error_std    = zeros(size(sigma_vals));
energy_error_mean = zeros(size(sigma_vals));
energy_error_std  = zeros(size(sigma_vals));

v_masked = v_noisy;  % masked signal already has dropout/quantization

n_realizations = 10;  % number of random trials per sigma
for k = 1:length(sigma_vals)
    sigma = sigma_vals(k);
    dist_err_realizations   = zeros(1, n_realizations);
    energy_err_realizations = zeros(1, n_realizations);
    
    for r = 1:n_realizations
        % Add Gaussian noise to velocity
        v_noisy_sigma = v_masked + sigma*randn(size(v_masked));
        v_noisy_sigma(isnan(v_masked)) = NaN;  % preserve masked points
        
        % Compute acceleration, power, energy, distance for noisy signal
        [accel_n, ~] = betterDifferentiator(t_noisy, v_noisy_sigma);
        P_n = eta .* min(Pclip, max(0, -m.*accel_n.*v_noisy_sigma));
        [E_n, ~] = betterIntegrator(t_noisy, P_n);
        [D_n, ~] = betterIntegrator(t_noisy, v_noisy_sigma);
        E_n_kWh = E_n / 3.6e6;
        
        % Compute percent error vs baseline
        dist_err_realizations(r)   = abs(D_n - D_base)/D_base*100;
        energy_err_realizations(r) = abs(E_n_kWh - E_base_kWh)/E_base_kWh*100;
    end
    
    % Store mean and standard deviation of error
    dist_error_mean(k)   = mean(dist_err_realizations);
    dist_error_std(k)    = std(dist_err_realizations);
    energy_error_mean(k) = mean(energy_err_realizations);
    energy_error_std(k)  = std(energy_err_realizations);
end

%% Determine max allowable sigma and output 
% Based on 2% distance error and 5% energy error thresholds
max_sigma_dist   = max(sigma_vals(dist_error_mean  <= 2));
max_sigma_energy = max(sigma_vals(energy_error_mean <= 5));
max_allowable_sigma = min(max_sigma_dist, max_sigma_energy);

% Find closest sigma index and associated standard deviation
[~, idx_sigma] = min(abs(sigma_vals - max_allowable_sigma));
sigma_uncertainty_dist   = dist_error_std(idx_sigma);
sigma_uncertainty_energy = energy_error_std(idx_sigma);
sigma_uncertainty = min(sigma_uncertainty_dist, sigma_uncertainty_energy);

% Print results
fprintf('Noise Analysis Results (mean ± std):\n');
fprintf('Max σ allowed (distance 2%%): %.3f m/s\n', max_sigma_dist);
fprintf('Max σ allowed (energy 5%%):   %.3f m/s\n', max_sigma_energy);
fprintf('Overall allowable σ: %.3f ± %.3f m/s\n', max_allowable_sigma, sigma_uncertainty);

%% Plots
%  Plot distance error vs sigma
figure;
fill([sigma_vals fliplr(sigma_vals)], ...
     [dist_error_mean+dist_error_std fliplr(dist_error_mean-dist_error_std)], ...
     [1 0.9 .9], 'EdgeColor','none'); hold on;
plot(sigma_vals, dist_error_mean, 'r','LineWidth',1.5); 
yline(2,'b--','2% Limit'); grid on;
xlabel('\sigma (m/s)'); ylabel('Distance Error (%)');
title('Distance Error vs Gaussian Noise (mean ± std)');

% plot energy error vs sigma
figure;
fill([sigma_vals fliplr(sigma_vals)], ...
     [energy_error_mean+energy_error_std fliplr(energy_error_mean-energy_error_std)], ...
     [1 0.9 0.9], 'EdgeColor','none'); hold on;
plot(sigma_vals, energy_error_mean, 'r','LineWidth',1.5); 
yline(5,'b--','5% Limit'); 
[~, idx_allowable] = min(abs(energy_error_mean - 5));
plot(sigma_vals(idx_allowable), energy_error_mean(idx_allowable), 'bo', 'MarkerFaceColor','b', 'MarkerSize',6);
text(sigma_vals(idx_allowable)+0.025, energy_error_mean(idx_allowable)+0.2, ...
    sprintf('\\sigma = %.3f m/s', sigma_vals(idx_allowable)), 'HorizontalAlignment','center');

grid on;
xlabel('\sigma (m/s)'); ylabel('Regen Energy Error (%)');
title('Energy Error vs Gaussian Noise (mean ± std)');



