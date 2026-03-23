%%
% Script: beamAnalysis_Chambers.m

% This script computs and visualizes the minimum beam thickness required to
% avoid a 75 Hz resosance frequency for a cantilever with a tip mass, along
% with visualizing the performance of three analysis methodologies. The
% characterist equation is based on a Euler–Bernoulli cantilever beam model
% which is evaluated across a sweep of thickness values, and solved for the
% nondimensional paramter of frequency, lambda, using the following root
% finding methods: Bisection, Newton-Raphson, and MATLAB's fzero. Lambda is
% then used as the final paramater needed to solve for the thickness given
% a rearranged frequency formula, with all other paramters needed given by
% the client. The smallest thickness that meets the target is identified.
% The method performance, based on speed, iteration count, and robustness
% from intial guesses is also evaluated and visualized. A further analysis
% is also done given a small range of payload variation, comparing how
% thickness needed changes based on the lumped tip mass changing. This is
% done to provide certainty in the stock aluminum choice given future
% design revisions.
% 
% Along with the raw computation of the models and visualations, the major 
% deliverables of this script are to provide a baseline for a recommended
% Aluminum 6061 stock, and also a comparison of analysis methodologies for
% the client to implement.

% The script calls the following functions:
% [x_r, iter, time] = bisectionMethod( f, x_l, x_u, e_s )
% [x_r, iter, time] = newtonRaphsonMethod( f, dfdx, x_0, e_s )
% e_s = significantDigitsPercentTolerance( n )

% Outputs: 
% - Figure 1: g(λ) with zero crossing and bounds
% - Figure 2: frequency vs thickness with target line and shaded regions
% - Figure 3: runtime vs specified percent tolerance (all methods)
% - Figure 4: iteration count vs specified percent tolerance (all methods)
% - Figure 5: bisection failure demonstration (non-bracketing example)
% - Figure 6: Newton–Raphson failure demonstration 1
% - Figure 7: Newton–Raphson failure demonstration 2
% - Figure 8 : fzero roots from different initial guesses
% - Figure 9: required thickness vs payload mass 
% - Figure 10: frequency response of a range of thickness for mass variation
% - Printed tables: initial guess vs root for Newton–Raphson and fzero

% Author:       Isaac Chambers
% Section:      ME 2016 - C
% Assignment:   MP2
% Date:         10/12/2025

clear, clc, close all; 

%% Part 1) Minimum Thickness Determination, 3 Methods used

% Initial Conditions 
L = .5;         % beam length (m)
b = .03;        % beam width (m)
E = 68.9e9;     % Young's Modulus (Pa)
rho = 2700;     % material density (kg/m^3)
M = 1.2;        % payload mass (kg)

thickness_test_range = linspace(0.02, 0.08, 1000);  % thickness sweep (m)
hz_target = 75;                                     % required natural frequency (Hz)

thickness_guess_lower_bound = .1;  % given from project instructions
thickness_guess_upper_bound = 3;   % given from project instructions
thickness_guess_midpoint = (thickness_guess_lower_bound + thickness_guess_upper_bound) / 2;

client_given_error_target = significantDigitsPercentTolerance(5);  % numeric tolerance from 5 sig figs, specified by client


% Bisection Method Analysis 
lambda_bisection = zeros(size(thickness_test_range));  % take care of initializations outside of loop for quicker runtime - consistent throughout
hz_bisection = zeros(size(thickness_test_range));

for i = 1:length(thickness_test_range)   % sweep thickness values, solve for lambda roots, convert to frequency, and find the minimum thickness that reaches the 75 Hz targets
    t = thickness_test_range(i);
    mu = M / (rho * b * t * L);     % payload-to-beam mass ratio per length

    % Characteristic equation for a cantilevered beam with tip mass (function of lambda)
    g = @(x) 1 + cosh(x).*cos(x) + mu.*x.*(sinh(x).*cos(x) - cosh(x).*sin(x));
   
    % Solve for lambda with bracketing
    [lambda_bisection(i), ~, ~] = bisectionMethod(g, thickness_guess_lower_bound, thickness_guess_upper_bound, client_given_error_target);

    I = b * t^3 / 12;               % second moment of area
    A = b * t;                      % cross-sectional area
    % Convert lambda to first natural frequency, rearranged from given equations
    hz_bisection(i) = (1/(2*pi)) * (lambda_bisection(i)^2 / L^2) * sqrt((E * I )/ (rho * A));
end

% First thickness meeting/exceeding the target frequency
idx_minimum_t = find(hz_bisection >= hz_target, 1, 'first');
t_min_b = thickness_test_range(idx_minimum_t);


% Newton Raphson Method Analysis
lambda_Newton_Raphson = zeros(size(thickness_test_range));
hz_Newton_Raphson = zeros(size(thickness_test_range));

for i = 1:length(thickness_test_range)
    t = thickness_test_range(i);
    mu = M / (rho * b * t * L);

    % Same characteristic equation and derivative
    g = @(x) 1 + cosh(x).*cos(x) + mu.*x.*(sinh(x).*cos(x) - cosh(x).*sin(x));
    dgdx = @(x) ((mu+1)*cos(x)-2*mu*x*sin(x))*sinh(x)+(-mu-1)*cosh(x)*sin(x); % verified with derivative-calculator.net
    
    % Solve for lambda with Newton–Raphson (requires derivative, uses midpoint of given bounds as a single initial guess)
    [lambda_Newton_Raphson(i), ~, ~] = newtonRaphsonMethod(g, dgdx, thickness_guess_midpoint, client_given_error_target);

    I = b * t^3 / 12;
    A = b * t;
    hz_Newton_Raphson(i) = (1/(2*pi)) * (lambda_Newton_Raphson(i)^2 / L^2) * sqrt((E * I )/ (rho * A));
end

% First thickness meeting/exceeding the target frequency
idx_minimum_t = find(hz_Newton_Raphson >= hz_target, 1, 'first');
t_min_NR = thickness_test_range(idx_minimum_t);


% Matlab fzero Function Analysis
lambda_fzero = zeros(size(thickness_test_range));
hz_fzero = zeros(size(thickness_test_range));

for i = 1:length(thickness_test_range)
    t = thickness_test_range(i);
    mu = M / (rho * b * t * L);

    g = @(x) 1 + cosh(x).*cos(x) + mu.*x.*(sinh(x).*cos(x) - cosh(x).*sin(x));
   
    % Solve for lambda with fzero (built in matlab root-finder function using a single initial guess)
    lambda_fzero(i) = fzero(g, thickness_guess_midpoint);
    
    I = b * t^3 / 12;
    A = b * t;
    hz_fzero(i) = (1/(2*pi)) * (lambda_fzero(i)^2 / L^2) * sqrt((E * I )/ (rho * A));
end

% First thickness meeting/exceeding the target frequency
idx_minimum_t = find(hz_fzero >= hz_target, 1, 'first');
t_min_fzero = thickness_test_range(idx_minimum_t);


%% Part 2) Performance Analysis 

t_minimum_value = t_min_fzero;   % establish minimum t value as evaluated in previous section; all values were the same, t_min_fzero used as fzero is a more advanced function

mu = M / (rho * b * t_minimum_value * L);
g = @(x) 1 + cosh(x).*cos(x) + mu.*x.*(sinh(x).*cos(x) - cosh(x).*sin(x));  % redefine g and dgdx outside of the loops 
dgdx = @(x) ((mu+1)*cos(x) - 2*mu*x*sin(x))*sinh(x) + (-mu-1)*cosh(x)*sin(x);  

number_of_sigfig_tolerance = linspace(2,10,1000);  % sweep tolerance over a range of desired significant digits (range and study resolution given) 

% initialize performance trackers for iteration count and runtime
bisection_iteration_count = zeros(size(number_of_sigfig_tolerance));
bisection_time = zeros(size(number_of_sigfig_tolerance));

Newton_Raphson_iteration_count = zeros(size(number_of_sigfig_tolerance));
Newton_Raphson_time = zeros(size(number_of_sigfig_tolerance));

fzero_iteration_count = zeros(size(number_of_sigfig_tolerance));
fzero_time = zeros(size(number_of_sigfig_tolerance));

error_values = zeros(size(number_of_sigfig_tolerance));

% Measure timing/iterations across tolerances
% Also, tic and toc are very volatile at the quickness of the runtimes for these functions (very fast), so the following takes medians over many repeitions to smooth out noise
for k = 1:length(number_of_sigfig_tolerance)

    current_error_target = significantDigitsPercentTolerance(number_of_sigfig_tolerance(k));
    error_values(k) = current_error_target;

    time_averaging_repetitions = 200;  % higher repetitions take longer to generate the figure, but provide a smoother look at the actual runtime trend of the different methods
    discard = 20;                      % trimmed first and last repetitions to attempt to account for/remove runtime spikes from the compilier 
    times_b = zeros(1, time_averaging_repetitions);
    times_NR = zeros(1, time_averaging_repetitions);
    times_f = zeros(1, time_averaging_repetitions);

    % Bisection average
    for r = 1:time_averaging_repetitions
        [~, bisection_iteration_count(k), this_time] = bisectionMethod(g, thickness_guess_lower_bound, thickness_guess_upper_bound, current_error_target);
        times_b(r) = this_time;
    end
    valid_b = times_b(discard+1 : end-discard);
    bisection_time(k) = median(valid_b);

    % Newton-Raphson average
    for r = 1:time_averaging_repetitions
        [~, Newton_Raphson_iteration_count(k), this_time] = newtonRaphsonMethod(g, dgdx, thickness_guess_midpoint, current_error_target);
        times_NR(r) = this_time;
    end
    valid_NR = times_NR(discard+1 : end-discard);
    Newton_Raphson_time(k) = median(valid_NR);

    % fzero average
    for r = 1:time_averaging_repetitions
        tic
        [~, ~, ~, output] = fzero(g, thickness_guess_midpoint);
        times_f(r) = toc;
    end
    valid_f = times_f(discard+1 : end-discard);
    fzero_time(k) = median(valid_f);
    fzero_iteration_count(k) = output.iterations;

end

% The following describes two tables showing the effect of a range of inital guesses on the root found by newton raphson and fzero respectively. These tables are used to compare robustness of the functions.

x0_guess_vals = linspace(0,10,50);   % arbitrary range of comparison
N = length(x0_guess_vals);

temp_root_Newton_Raphson = zeros(N,1);
temp_root_fzero = zeros(N,1);

for i = 1:N

    % Newton–Raphson 
    [x_r, ~, ~] =newtonRaphsonMethod(g, dgdx, x0_guess_vals(i), current_error_target);
    temp_root_Newton_Raphson(i) = x_r;

    % fzero
    [temp_root_fzero(i)] = fzero(g, x0_guess_vals(i));

end

table_Newton_Raphson = table(x0_guess_vals(:), temp_root_Newton_Raphson, ...
    'VariableNames', {'x0','root_NR'});

table_fzero = table(x0_guess_vals(:), temp_root_fzero, ...
    'VariableNames', {'x0','root_fzero'});

disp('--- Newton–Raphson results ---');
disp(table_Newton_Raphson);

disp('--- fzero results ---');
disp(table_fzero);

%% Figures for presentation 

% g(x) zoomed in on where lamda crosses x axis given t_min
x_vals = linspace(0, thickness_guess_upper_bound, 1000);   
g_vals = g(x_vals);

figure;
plot(x_vals, g_vals, 'LineWidth', 1.5);
hold on;
yline(0, 'k--', 'LineWidth', 1);  
plot(thickness_guess_lower_bound, g(thickness_guess_lower_bound), 'ro', 'MarkerFaceColor', 'r');
plot(thickness_guess_upper_bound, g(thickness_guess_upper_bound), 'go', 'MarkerFaceColor', 'g');
plot(x_vals(abs(g_vals)==min(abs(g_vals))), 0, 'ks', 'MarkerFaceColor', 'k');
xlabel('\lambda');
ylabel('g(\lambda)');
title('Characteristic Equation g(\lambda) vs \lambda');
grid on;


% Thickness vs frequency
figure; 
hold on; 
plot(thickness_test_range*1000, hz_Newton_Raphson, 'LineWidth', 1.5); % *1000 to plot in mm
area1 = linspace(min(thickness_test_range),t_minimum_value,length(hz_Newton_Raphson(1:idx_minimum_t))).*1000;
area(area1,hz_Newton_Raphson(1:idx_minimum_t), 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none')
area2 = linspace(t_minimum_value,max(thickness_test_range), length(hz_Newton_Raphson(idx_minimum_t:end))).*1000;
area(area2,hz_Newton_Raphson(idx_minimum_t:end), 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none')
tlabel = sprintf('Minimum Thickness Required\n t ≥ %2.3f mm', t_minimum_value*1000);
yline(hz_target, 'k--', 'Freqeuncy Target, 75 Hz', ...
      'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'middle', ...    ]
      'LabelOrientation', 'horizontal');

plot(thickness_test_range(idx_minimum_t)*1000, hz_Newton_Raphson(idx_minimum_t), 'ko', 'MarkerFaceColor', 'k'); % *1000 to plot in mm
text(thickness_test_range(idx_minimum_t)*1000+7, hz_Newton_Raphson(idx_minimum_t)-10, tlabel, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'Rotation', 0);
xlabel('Thickness t (mm)'); ylabel('First Natural Frequency f_1 (Hz)');
title('Beam First Mode Frequency vs. Thickness');
grid on;


% Runtime vs tolerance
figure;
plot(error_values, bisection_time, '-r', ...
     error_values, Newton_Raphson_time, '-g', ...
     error_values, fzero_time, '-b', 'LineWidth', 1.2);
set(gca, 'XDir', 'reverse');
xlabel('Specified Percent Tolerance (%)');
ylabel('Runtime (seconds)');
legend('Bisection','Newton–Raphson','fzero','Location','best');
title('Runtime vs Specified Percent Tolerance');
grid on;

% Iterations count vs tolerance
figure;
plot(error_values, bisection_iteration_count, '-r', ...
     error_values, Newton_Raphson_iteration_count, '-g', ...
     error_values, fzero_iteration_count, '-b', 'LineWidth', 1.2);
set(gca, 'XDir', 'reverse');
xlabel('Specified Percent Tolerance (%)');
ylabel('Iteration Count');
legend('Bisection','Newton–Raphson','fzero','Location','best');
title('Iteration Count vs Specified Percent Tolerance');
grid on;

% Bisection failure case
x_l_bad = 2.5; x_u_bad = 3.0; % interval not containing first root
x_vals = linspace(0, 4, 1000);
figure;
plot(x_vals, g(x_vals), 'LineWidth', 1.5); hold on;
yline(0, 'k--');
xline(x_l_bad, 'r--', 'Lower bound');
xline(x_u_bad, 'g--', 'Upper bound');
title('Bisection: Improper Bounds Fail to Bracket Root');
xlabel('\lambda'); ylabel('g(\lambda)');
legend('g(\lambda)','Zero line','x_l','x_u');
grid on;

% Newton failure case
x_bad_newton_guess = 3.5;  
max_iter = 10;
iterVals = zeros(1, max_iter);
iterVals(1) = x_bad_newton_guess;

for k = 2:max_iter
    iterVals(k) = iterVals(k-1) - g(iterVals(k-1)) / dgdx(iterVals(k-1));
end

x_vals = linspace(0, 7.5, 500);
figure;
plot(x_vals, g(x_vals), 'b-', 'LineWidth', 1.5); hold on;
yline(0, 'k--');
plot(iterVals, g(iterVals), 'ro-', 'MarkerFaceColor', 'r', 'DisplayName', 'Iterations');
xlabel('\lambda'); ylabel('g(\lambda)');
title('Newton–Raphson Iteration Path (Divergence Example)');
legend('g(\lambda)', 'Zero line', 'Iterations');
grid on;

text(iterVals(1), g(iterVals(1)), '  Start', 'Color', 'r', 'FontWeight', 'bold');
text(iterVals(end), g(iterVals(end)), '  After 10 iters', 'Color', 'r', 'FontWeight', 'bold');

figure;
plot(1:max_iter, iterVals, 'ro-', 'LineWidth', 1.5);
xlabel('Iteration Number');
ylabel('\lambda estimate');
title('Newton–Raphson Iteration Sequence (Divergence Example)');
grid on;


%% Fzero failure case
x_guesses = [-3, thickness_guess_midpoint, 10];
colors = ['r','g','b'];
figure; hold on;
fplot(g, [-6 10.3], 'k-', 'LineWidth', 1.5);

for i = 1:length(x_guesses)
    root = fzero(g, x_guesses(i));
    plot(root, 0, 'o', 'MarkerFaceColor', colors(i), 'DisplayName', sprintf('Start %.1f -> Root %.3f', x_guesses(i), root));
end

xlabel('\lambda'); ylabel('g(\lambda)');
title('fzero Root Selection Dependence on Initial Guess');
legend show; grid on;





%% Further Analysis on the stock recomendation for future revision robustness

M_span = linspace(0.7*M, 1.5*M, 25); % mass variation range
t_req = zeros(size(M_span));

for j = 1:length(M_span)
    Mj = M_span(j);
    found = false;
    for i = 1:length(thickness_test_range)  % t_range in meters
        t = thickness_test_range(i);
        I = b * t^3 / 12;
        A = b * t;
        mu = Mj / (rho * b * t * L);
        g_mass_var = @(x) 1 + cosh(x).*cos(x) + mu.*x.*(sinh(x).*cos(x) - cosh(x).*sin(x));
        lam = fzero(g_mass_var, thickness_guess_midpoint);   
        hz = (1/(2*pi)) * (lam^2 / L^2) * sqrt((E * I )/ (rho * A));
        if hz >= hz_target
            t_req(j) = t;            % minimal t meeting 75 Hz for this M
            found = true;
            break
        end
    end
end

t_max_stock = 50.8e-3;  % 50.8 mm
figure;
plot(M_span, t_req*1000, 'b-', 'LineWidth', 1.6); % *1000 to plot in mm
hold on;
yline(t_max_stock*1000, 'k--', 't_{stock} = 50.8 mm', 'LabelHorizontalAlignment','left');
xlabel('Payload Mass M (kg)');
ylabel('Required Thickness t_{req} (mm) to meet 75 Hz');
title('Required Thickness vs Payload Mass');
grid on;


% Mass variation analysis: +/-10% tip mass
M_nom = M;                  % original mass (1.20 kg)
M_low = 0.9 * M_nom;
M_high = 1.1 * M_nom;

t_vals = thickness_test_range;           % same thickness range
hz_low = zeros(size(t_vals));
hz_nom = zeros(size(t_vals));
hz_high = zeros(size(t_vals));

for i = 1:length(t_vals)
    t = t_vals(i);
    
    I = b * t^3 / 12;
    A = b * t;
    
    % Lambda solving for each mass level using fzero
    % g depends on mu => different for each mass case

    % Low mass
    mu_low = M_low / (rho * b * t * L);
    g_low = @(x) 1 + cosh(x).*cos(x) + mu_low.*x.*(sinh(x).*cos(x) - cosh(x).*sin(x));
    lam_low = fzero(g_low, thickness_guess_midpoint);
    hz_low(i) = (1/(2*pi)) * (lam_low^2 / L^2) * sqrt((E * I )/ (rho * A));
    
    % Nominal mass
    mu_nom = M_nom / (rho * b * t * L);
    g_nom = @(x) 1 + cosh(x).*cos(x) + mu_nom.*x.*(sinh(x).*cos(x) - cosh(x).*sin(x));
    lam_nom = fzero(g_nom, thickness_guess_midpoint);
    hz_nom(i) = (1/(2*pi)) * (lam_nom^2 / L^2) * sqrt((E * I )/ (rho * A));

    % High mass
    mu_high = M_high / (rho * b * t * L);
    g_high = @(x) 1 + cosh(x).*cos(x) + mu_high.*x.*(sinh(x).*cos(x) - cosh(x).*sin(x));
    lam_high = fzero(g_high, thickness_guess_midpoint);
    hz_high(i) = (1/(2*pi)) * (lam_high^2 / L^2) * sqrt((E * I )/ (rho * A));
end


figure;
plot(t_vals*1000, hz_nom, 'b-', 'LineWidth', 1.6); hold on;
plot(t_vals*1000, hz_low, 'g--', 'LineWidth', 1.2);
plot(t_vals*1000, hz_high, 'r--', 'LineWidth', 1.2);
yline(75, 'k--', '75 Hz Requirement', 'LabelHorizontalAlignment','left');

xlabel('Thickness t (mm)');
ylabel('First Natural Frequency f_1 (Hz)');
title(' Frequency Response of ±10% Mass Variation');
legend('Nominal Mass','-10% Mass', '+10% Mass');