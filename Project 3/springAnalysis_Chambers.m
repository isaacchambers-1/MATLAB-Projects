%%
% Script: springModelComparison_Chambers.m
%
% This script performs regression analysis and model validation for a 
% nonlinear spring system tested under two experimental setups: a drive 
% test and a lab test. The objective is to identify the most appropriate 
% model for characterizing the force–displacement relationship, evaluate 
% model accuracy, and visualize both interpolation and extrapolation 
% behaviors.
%
% The models compared include linear, quadratic, cubic, power, logarithmic, 
% and exponential regressions, fitted using a least-squares regression 
% function. Each model's residual distribution and root-mean-square error 
% (RMSE) are computed and compared to determine the best overall fit for 
% the nonlinear spring data. A cubic model is selected for further analysis 
% and extrapolation.
%
% In the second section, cubic spline interpolation is applied to the lab 
% data and compared with the cubic regression model using both absolute 
% and relative error metrics. The third section extends both models beyond 
% the measured displacement range to assess predictive stability and 
% divergence in the extrapolated region.
%
% The script calls the following functions:
% [y_fit, residuals] = lsRegression(x_data, y_data, fitType)
% rmse = rmsError(y_true, y_pred)
%
% Outputs:
% - Figure 1: Lab test Force vs. Displacement
% - Figure 2: Drive test Force vs. Displacement
% - Figure 3: All regression model fits vs. drive test data
% - Figure 4: Selected regression model fits (subset)
% - Figure 5: All regression residuals
% - Figure 6: Selected regression residuals
% - Figure 7: RMSE comparison bar chart
% - Figure 8: Spline interpolation of lab data
% - Figure 9: Comparison between spline and cubic regression
% - Figure 10: Spline vs regression with error overlays
% - Figure 11: Relative error between models
% - Figure 12: Extrapolation of models up to 0.4 m displacement
%
% Author:       Isaac Chambers
% Section:      ME 2016 - C
% Assignment:   MP3
% Date:         10/29/2025


clear, clc, close all;

load("drivetest.mat")
load("labtest.mat")


%% Part 1 

figure;
scatter(x_lab, F_lab, 'r', 'filled')
xlabel('Displacement (m)') 
ylabel('Force (N)')
title('Lab test Force vs. Displacement')

figure;
scatter(x_drive, F_drive, 'r', 'filled')
xlabel('Displacement (m)') 
ylabel('Force (N)')
title('Drive Test Force vs. Displacement')

figure
scatter(x_drive, F_drive, 25, 'r', 'filled')
xlabel('Displacement (m)')
ylabel('Force (N)')
hold on
grid on

fit_types = {'linear','quadratic','cubic','power', 'log', 'exponential'};

for i = 1:numel(fit_types)
    [y_line, ~] = lsRegression(x_drive, F_drive, fit_types{i});
    plot(x_drive, y_line, 'LineWidth', 1.5, 'DisplayName', fit_types{i})
end

legend('Data', fit_types{:}, 'Location', 'best')
title('All Regression Model Fits vs. Drive Test Data')

figure
scatter(x_drive, F_drive, 25, 'r', 'filled')
xlabel('Displacement (m)')
ylabel('Force (N)')
hold on
grid on

x_line = linspace(min(x_drive), max(x_drive), 500)';
fit_types_select = {'linear','quadratic','cubic','power'};

for i = 1:numel(fit_types_select)
    [y_line, ~] = lsRegression(x_drive, F_drive, fit_types_select{i});
    plot(x_drive, y_line, 'LineWidth', 2, 'DisplayName', fit_types_select{i})
end

legend('Data', fit_types_select{:}, 'Location', 'best')
title('Select Regression Model Fits vs. Drive Test Data')

figure
hold on
grid on

for i = 1:numel(fit_types)
    [~, residuals] = lsRegression(x_drive, F_drive, fit_types{i});
    scatter(x_drive, residuals, 'filled')
end
yline(0, 'LineWidth', 2.5, 'Color', 'black')
legend(fit_types{:})
title('All Regression Model Residuals vs. Displacement')

figure
hold on
grid on

for i = 1:numel(fit_types_select)
    [~, residuals] = lsRegression(x_drive, F_drive, fit_types_select{i});
    scatter(x_drive, residuals, 'filled')
end
yline(0, 'LineWidth', 2.5, 'Color', 'black')
legend(fit_types_select{:})
title('Select Regression Model Residuals vs. Displacement')

rmse_vals = zeros(numel(fit_types),1);
for i = 1:numel(fit_types)
    [y,resid] = lsRegression(x_drive, F_drive, fit_types{i});
    rmse_vals(i) = rmsError(F_drive, y);
end

T = table(char(fit_types'), rmse_vals, 'VariableNames', {'Model','RMSE'});
disp(T)

figure
bar(fit_types, rmse_vals)
ylabel('RMSE')
title('Goodness-of-Fit Comparison for All Models')
grid on

% cubic equation
[y_fit_best, ~] = lsRegression(x_drive, F_drive, 'cubic');
x = x_drive(:);
A = [x.^3 x.^2 x ones(size(x))];
b = pinv(A) * F_drive(:);
fprintf('Cubic model: F = %.4f*x^3 + %.4f*x^2 + %.4f*x + %.4f\n', b)


%% Part 2 
x_fine = linspace(min(x_lab), max(x_lab), 500)';

S = spline(x_lab, F_lab);
F_spline = ppval(S, x_fine);

[y_reg, ~] = lsRegression(x_drive, F_drive, 'cubic');
F_reg_interp = interp1(x_drive, y_reg, x_fine, 'linear','extrap');


figure
scatter(x_lab, F_lab, 40, 'r', 'filled')
hold on
plot(x_fine, F_spline, 'b-', 'LineWidth', 2)
xlabel('Displacement (m)')
ylabel('Force (N)')
title('Cubic Spline Interpolation of Lab Test Data')
legend('Lab Data','Cubic Spline','Location','best')
grid on

figure
plot(x_fine, F_spline, 'b','LineWidth',2)
hold on
plot(x_fine, F_reg_interp, 'r--','LineWidth',2)
scatter(x_lab, F_lab,40,'k','filled')
xlabel('Displacement (m)')
ylabel('Force (N)')
legend('Spline Interpolation','Cubic Regression','Lab Data','Location','best')
title('Comparison: Spline vs Cubic Regression')
grid on


error_abs = abs(F_spline - F_reg_interp);
error_rel = 100 * abs((F_spline - F_reg_interp) ./ F_spline);

figure
subplot(2,1,1)
yyaxis left
plot(x_fine, F_spline, 'b-', 'LineWidth', 2)
hold on
plot(x_fine, F_reg_interp, 'r--', 'LineWidth', 2)
scatter(x_lab, F_lab, 40, 'k', 'filled')
ylabel('Force (N)')
yyaxis right
plot(x_fine, error_abs, 'm-', 'LineWidth', 1.8)
ylabel('|Error| (N)')
xlabel('Displacement (m)')
title('Spline vs Regression with Absolute Error Overlay')
legend('Spline','Regression','Lab Data','|Error|','Location','best')
grid on

subplot(2,1,2)
plot(x_fine, error_rel, 'c-', 'LineWidth', 2)
xlabel('Displacement (m)')
ylabel('Relative Error (%)')
title('Relative Error Between Models')
ylim([0, 50])
grid on


%% Part 3 

x_ext = linspace(0, 0.4, 500)';

S = spline(x_lab, F_lab);
F_spline_ext = ppval(S, x_ext);

[y_reg, ~] = lsRegression(x_drive, F_drive, 'cubic');
F_reg_ext = interp1(x_drive, y_reg, x_ext, 'linear', 'extrap');

figure
hold on
scatter(x_lab, F_lab, 40, 'k', 'filled')
plot(x_ext, F_spline_ext, 'b-', 'LineWidth', 2)
plot(x_ext, F_reg_ext, 'r--', 'LineWidth', 2)
xlabel('Displacement (m)')
ylabel('Force (N)')
title('Model Extrapolation to 0.4 m')
legend('Lab Data','Spline Interpolation','Cubic Regression','Location','best')
grid on



figure
plot(x_ext, F_reg_ext, 'r--', 'LineWidth', 2)
hold on
plot(x_ext, F_spline_ext, 'b-', 'LineWidth', 1.5)
xline(0.3, 'k--', 'LineWidth', 1)
xlabel('Displacement (m)')
ylabel('Force (N)')
title('Extrapolation Beyond Measured Range (0.3 m to 0.4 m)')
legend('Cubic Regression','Spline Interpolation', 'Data Limit','Location','best')
grid on
xlim([0.25 0.4])
