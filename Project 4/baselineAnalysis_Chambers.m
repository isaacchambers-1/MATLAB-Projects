%% Script: noiseBaselineAnalysis_Chambers.m
%
% This script performs a baseline analysis of regenerative braking
% power recovery using UDDS velocity data. The workflow includes:
%   - Loading the UDDS dataset
%   - Converting velocity units from mph to m/s
%   - Computing acceleration using a numerical differentiator
%   - Calculating instantaneous regen power based on vehicle mass,
%     efficiency, and a power clip limit
%   - Computing cumulative distance traveled and total regen energy
%   - Visualizing velocity, acceleration, regen power, and cumulative energy
% The overall goal is to provide an ideal baseline analysis given trusted
% data. 
%
% The script also prints method usage counts for integrators and differentiator.
%
%The script calls the following functions:
% [D, counts] = betterDifferentiator(t, v)
%   Computes numerical derivative of a velocity signal to get acceleration.
%   Returns derivative vector and counts of finite difference methods used.
%
% [I, counts] = betterIntegrator(t, x)
%   Computes numerical integral of a signal over time.
%   Returns integrated value and counts of integration methods used.
%
% Outputs:
%   - Figures:
%       1) UDDS velocity profile
%       2) Acceleration vs time
%       3) Instantaneous regen power vs time
%       4) Cumulative regen energy vs time
%   - Console output:
%       Total distance traveled, total regen energy (kWh)
%       Integrator and differentiator method counts
%
% Note: This script also loads 'udds_LA4_data' from 'uddscol.mat',
%       which must be in the MATLAB path or current folder.
%
% Author:       Isaac Chambers
% Section:      ME 2016 - C
% Assignment:   MP4
% Date:         11/23/2025

clc;        
format longG;               
%% Load and unit convert UDDS data
load uddscol              
udds_time = udds_LA4_data(:,1);            
udds_vel  = (1609.344/3600) * udds_LA4_data(:,2); 

%% Computations
m      = 1500;             % Vehicle mass (kg)
eta    = 0.75;             % Regenerative braking efficiency
P_clip = 80000;            % Maximum regen power allowed (W)

[total_distance, dist_counts] = betterIntegrator(udds_time, udds_vel); % integrate velocity over time
[udds_accel, ~] = betterDifferentiator(udds_time, udds_vel); % acceleration from velocity

P_regen = (eta * min(P_clip, max(0, -m .* udds_accel .* udds_vel))) * (1/1000); % kW

[Energy_regen, regen_counts] = betterIntegrator(udds_time, P_regen); % integrate power to get total energy
Energy_regen_kWh = Energy_regen / 3600;    % convert kWs to kWh

%% Display results
disp("Baseline Results:")
fprintf("Total Distance Traveled: %.3f km\n", total_distance/1000);
fprintf("Total Regen Energy (kWh): %.4f kWh\n", Energy_regen_kWh);

disp("Method Usage Counts:")
disp("Distance integrator counts [trap, 1/3, 3/8]:")
disp(dist_counts')

disp("Regen energy integrator counts [trap, 1/3, 3/8]:")
disp(regen_counts')

disp("Differentiator counts [CDD, Interp-DD]:")
[~, diff_counts] = betterDifferentiator(udds_time, udds_vel);
disp(diff_counts')


%% Plots
% Plot velocity & accel vs time
figure;
yyaxis left
plot(udds_time, udds_vel, 'b', 'LineWidth', 1.5);
ylabel('Velocity [m/s]');
xlabel('Time [s]');
grid on;

yyaxis right
plot(udds_time, udds_accel, 'r', 'LineWidth', 1.4);
ylabel('Acceleration [m/s^2]');

title('UDDS Velocity and Acceleration (betterDifferentiator)');
legend('Velocity','Acceleration','Location','best');


%Differentiaor Showcase
x = [0 1 2 3.2 4.5 5.5 6 7.1 8];
y = sin(x);

[D, counts] = betterDifferentiator(x, y); 

uniform_idx = [true true true false false true true false true];
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];

figure('Color','w'); hold on; grid on;

xx = linspace(min(x), max(x), 200);
yy = sin(xx);
plot(xx, yy, 'k--','LineWidth',1.5,'DisplayName','Original Function');

plot(x, y, 'ko','MarkerFaceColor','k','DisplayName','Sample Points');

for i = 1:length(x)-2
    if uniform_idx(i)
        plot(x(i:i+2), y(i:i+2), '-o','Color',colors(1,:),'LineWidth',2, 'MarkerFaceColor',colors(1,:));
        text(mean(x(i:i+2)), mean(y(i:i+2))+0.1, 'CDD','Color',colors(1,:),'FontWeight','bold','HorizontalAlignment','center');
    else
        plot(x(i:i+2), y(i:i+2), '-o','Color',colors(2,:),'LineWidth',2, 'MarkerFaceColor',colors(2,:));
        text(mean(x(i:i+2)), mean(y(i:i+2))-0.15, 'IPDD','Color',colors(2,:),'FontWeight','bold','HorizontalAlignment','center');
    end
end

xlabel('x'); ylabel('y'); 
title('betterDifferentiator: CDD vs IPDD');
legend('Original Function: sin(x)','Sample Points','Location','best');
ylim([-1.5 1.5]);


% Results Plot
figure;
hold on; 
grid on;

yyaxis left
plot(udds_time, P_regen, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Instantaneous Regen Power [kW]');
ylabel('Power [kW]');

yyaxis right
cumulative_energy = cumsum(P_regen .* [0; diff(udds_time)]); % integrate power to get cumulative energy
plot(udds_time, cumulative_energy/3600, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Cumulative Energy [kWh]');
ylabel('Energy [kWh]');

xlabel('Time [s]');
title('Baseline Regenerative Braking: Power & Cumulative Energy');
legend('Location','best');



%% Validation of betterIntegrator vs trapz
format short

% Define test functions and analytical integrals
funcs = {@(x) x.^2, @(x) sin(x), @(x) exp(-x)};
exact_integrals = [1/3, 1-cos(1), 1-exp(-1)]; % integral from 0 to 1

x = linspace(0,1,11); % uniform spacing
results = zeros(length(funcs),4); % [betterIntegrator, trapz, %error_better, %error_trapz]

for k = 1:length(funcs)
    y = funcs{k}(x);
    
    % BetterIntegrator
    [I_b, counts] = betterIntegrator(x, y);
    
    % trapz
    I_t = trapz(x,y);
    
    % Store results
    results(k,1) = I_b;
    results(k,2) = I_t;
    results(k,3) = abs(I_b - exact_integrals(k))/exact_integrals(k) * 100;
    results(k,4) = abs(I_t - exact_integrals(k))/exact_integrals(k) * 100;
end

% Display table
T = table({'x^2';'sin(x)';'exp(-x)'}, results(:,1), results(:,2), results(:,3), results(:,4), ...
    'VariableNames',{'Function','BetterIntegrator','Trapz','Error_Better_%','Error_Trapz_%'});
disp(T);