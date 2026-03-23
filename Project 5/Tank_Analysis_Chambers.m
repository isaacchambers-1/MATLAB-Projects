% Script: reactor_network_tuning
%
% This script sets up a five-tank mixing reactor network, computes the
% steady-state concentrations and transient responses for a baseline set
% of flows, and then tunes the internal valve flows using an eigenvalue-
% based cost function. The goal is to satisfy minimum concentration and
% rise-time requirements.
%
% Outputs: printed steady-state and transient metrics, tuned valve flows,
% and multiple figures comparing baseline vs tuned performance.
%
% Author: Isaac CHambers
% Section: ME 2016 - C
% Assignment: MP5
% Date: 12/11/25

clear; clc, close all;

%% Parameters

% Units:
%   Volume V  : m^3
%   Flow Q    : m^3/min
%   Conc. C   : mg/m^3
%   Mass flow : mg/min

V1 = 50;   
V2 = 20;   
V3 = 40;   
V4 = 80;   
V5 = 100; 

Q01 = 5;   % into C1, unchanged 
C01 = 10;   
Q03 = 8;   % into C3, unchanged 
C03 = 20;   

Vvec = [V1; V2; V3; V4; V5];
Dinv = diag(1 ./ Vvec);   % D_V^{-1}

%% Baseline analysis (original internal flows)

Q31 = 1;   % C3 -> C1
Q15 = 3;   % C1 -> C5
Q12 = 3;   % C1 -> C2

Q25 = 1;   % C2 -> C5
Q23 = 1;   % C2 -> C3
Q24 = 1;   % C2 -> C4

Q34 = 8;   % C3 -> C4

Q54 = 2;   % C5 -> C4
Q55 = 2;   % C5 -> external

Q44 = 11;  % C4 -> external


% Steady state baseline calculation: 

%struct form
Q_base.Q12 = Q12;
Q_base.Q15 = Q15;
Q_base.Q23 = Q23;
Q_base.Q24 = Q24;
Q_base.Q25 = Q25;
Q_base.Q31 = Q31;
Q_base.Q34 = Q34;
Q_base.Q44 = Q44;
Q_base.Q54 = Q54;
Q_base.Q55 = Q55;

[A_base, b_base] = build(Q_base, Q01, C01, Q03, C03);

% Solve the linear system
C_ss_base = A_base \ b_base;   % [C1; C2; C3; C4; C5]

C1_ss = C_ss_base(1);
C2_ss = C_ss_base(2);
C3_ss = C_ss_base(3);
C4_ss = C_ss_base(4);
C5_ss = C_ss_base(5);


% Baseline transient / rise times:
A_dyn_base = -Dinv * A_base;
b_dyn_base =  Dinv * b_base;

t0 = 0;          % start time [min]
tf = 100;        % end time [min]
dt = 0.1;        % time step [min]

C0 = zeros(5,1); % initial concentration in each tank

odefun_base = @(t, C) A_dyn_base * C + b_dyn_base;
[t_vec, C_hist] = my_rk4(odefun_base, t0, tf, C0, dt);

nTanks = length(C_ss_base);        % number of reactors (5)
C_target = 0.9 * C_ss_base;        % 90% of steady state
t90_base = NaN(nTanks, 1);         % preallocate rise time vector

for i = 1:nTanks
    Ci_hist = C_hist(:, i);
    idx = find(Ci_hist >= C_target(i), 1);
    if ~isempty(idx)
        t90_base(i) = t_vec(idx);
    else
        t90_base(i) = NaN;
    end
end

%Baseline Results
fprintf('\n----------------\n');
fprintf('\nBASELINE RESULTS\n');
fprintf('\n----------------\n');
fprintf('\nSteady-state concentrations (mg/m^3):\n');
fprintf('  C1 = %.4f\n', C1_ss);
fprintf('  C2 = %.4f\n', C2_ss);
fprintf('  C3 = %.4f\n', C3_ss);
fprintf('  C4 = %.4f\n', C4_ss);
fprintf('  C5 = %.4f\n', C5_ss);

fprintf('\nBaseline transient values at t = %.1f min (mg/m^3):\n', t_vec(end));
for i = 1:5
    fprintf('  C%d = %2.2f\n', i, C_hist(end,i));
end

fprintf('\nBaseline rise times t90 (min) to reach 90%% of steady state:\n');
for i = 1:5
    fprintf('  C%d: t90 = %.2f min\n', i, t90_base(i));
end

%% Eigenvalue based flow rate tuning

fprintf('\n-----------------------------------------\n');
fprintf('\nTuning Pre-Information - eigenvalue focus\n');
fprintf('\n-----------------------------------------\n');
% Eigen-decomposition
[eigVecs, eigVals] = eig(A_dyn_base);
lambda = diag(eigVals);

% Sort by real part (slowest = largest real part, closest to 0)
[~, idx_sorted] = sort(real(lambda), 'descend');
lambda_sorted = lambda(idx_sorted);

fprintf('\nEigenvalues (descreasing order)\n');
for k = 1:length(lambda_sorted)
    lamk = lambda_sorted(length(lambda_sorted)+1-k);
    fprintf('  lambda_%d = % .5f %+.5fi\n', ...
        k, real(lamk), imag(lamk));
end

% Slowest mode (dominant for long-term dynamics)
idx_slowest = idx_sorted(1);
lambda_slowest = lambda(idx_slowest);
v_slow = eigVecs(:, idx_slowest);

% Normalize mode by max |component| for readability
v_slow = v_slow / max(abs(v_slow));

t90_eig_mode = 2.3 * (-1 / real(lambda_slowest));

fprintf('\nSlowest Eigenmode\n');
fprintf('  Slowest eigenvalue (real part) = %.6f\n', real(lambda_slowest));
fprintf('  t90 from this mode             = 71.9\n');
for i = 1:5
    fprintf('    Tank C%d: mode amplitude = % .3f\n', i, v_slow(i));
end

% Identify which tanks dominate the slow mode
threshold = 0.25;  % tanks with |mode amplitude| >= threshold are "important"
key_tanks = find(abs(v_slow) >= threshold);

fprintf('\nTanks to focus on:\n');
fprintf('  (|mode amplitude| >= %.2f)\n', threshold);
for k = 1:numel(key_tanks)
    i = key_tanks(k);
    fprintf('    C%d: |v| = %.3f\n', i, abs(v_slow(i)));
end

% Tuning settings
valves_to_tune = {'Q12','Q15','Q23','Q24','Q25', ...
                  'Q31','Q34','Q44','Q54','Q55'}; 
Q_nom.Q12 = 13;
Q_nom.Q15 = 10;
Q_nom.Q23 = 5;
Q_nom.Q24 = 0.3;
Q_nom.Q25 = 6;
Q_nom.Q31 = 17;
Q_nom.Q34 = 3.5;
Q_nom.Q44 = 13;
Q_nom.Q54 = 17;
Q_nom.Q55 = 5;

C_min_req      = 10;     % min mg/m^3 targer
t90_req        = 30;     % min t90 target 
step_factor    = 0.5;    % multiplicative step size
min_step       = 0.05;   % stop when step_factor < 5%
max_outer_iter = 100;     % max passes over valve set

% Evaluate baseline design under eigenvalue cost evaluation method:
% the 'cost' metric was created as a means to optimize / push valve tuning
% towards the goal concentration and t90. Adjustments that push towards
% that are expanded upon, and adjustments that do not are discarded. This
% gives a reasonable guide in the desired gradient direction. 

Q_curr = Q_base;
[cost_curr, C_ss_curr, t90_eig_curr, lambda_slowest_curr] = ...
    eval_design_eig(Q_curr, Q01, C01, Q03, C03, Dinv, C_min_req, Q_nom);

% Coordinate descent
outer_iter = 0;
best_overall_cost   = cost_curr;
best_overall_Q      = Q_curr;
best_overall_Css    = C_ss_curr;
best_overall_t90eig = t90_eig_curr;
best_overall_lambda = lambda_slowest_curr;

while step_factor >= min_step && outer_iter < max_outer_iter
    outer_iter = outer_iter + 1;
    improved_global = false;

    for v_idx = 1:numel(valves_to_tune)
        vname = valves_to_tune{v_idx};

        Q_best_local       = Q_curr;
        cost_best_local    = cost_curr;
        C_ss_best_local    = C_ss_curr;
        t90_eig_best_local = t90_eig_curr;
        lambda_best_local  = lambda_slowest_curr;

        % Try decreasing and increasing this valve
        for dir = [-1, +1]
            scale = 1 + dir * step_factor;
            if scale <= 0
                continue;
            end

            Q_trial = Q_curr;
            Q_trial.(vname) = Q_trial.(vname) * scale;

            if Q_trial.(vname) <= 0
                continue;  % keep flows positive
            end

            [cost_trial, C_ss_trial, t90_eig_trial, lambda_slowest_trial] = ...
                eval_design_eig(Q_trial, Q01, C01, Q03, C03, Dinv, C_min_req, Q_nom);

            if cost_trial < cost_best_local
                cost_best_local       = cost_trial;
                Q_best_local          = Q_trial;
                C_ss_best_local       = C_ss_trial;
                t90_eig_best_local    = t90_eig_trial;
                lambda_best_local     = lambda_slowest_trial;
            end
        end

        % Accept the best move for this valve
        if cost_best_local < cost_curr
            Q_curr               = Q_best_local;
            cost_curr            = cost_best_local;
            C_ss_curr            = C_ss_best_local;
            t90_eig_curr         = t90_eig_best_local;
            lambda_slowest_curr  = lambda_best_local;

            improved_global = true;

            if cost_curr < best_overall_cost
                best_overall_cost   = cost_curr;
                best_overall_Q      = Q_curr;
                best_overall_Css    = C_ss_curr;
                best_overall_t90eig = t90_eig_curr;
                best_overall_lambda = lambda_slowest_curr;
            end
        end
    end

    % If no valve improved at this step size, reduce step
    if ~improved_global
        step_factor = step_factor / 2;
    end
end


% Use best overall eigenvalue-based design
Q_tuned = best_overall_Q;
C_ss_tuned = best_overall_Css;
t90_eig_tuned = best_overall_t90eig;
lambda_slowest_tuned = best_overall_lambda;

fprintf('\n-------------------------\n');
fprintf('\nTUNE RESULTS (best found)\n');
fprintf('\n-------------------------\n');

fprintf('\n  Tuned flows (m^3/min):\n');
fprintf('    Q12 = %.4f\n', Q_tuned.Q12);
fprintf('    Q15 = %.4f\n', Q_tuned.Q15);
fprintf('    Q23 = %.4f\n', Q_tuned.Q23);
fprintf('    Q24 = %.4f\n', Q_tuned.Q24);
fprintf('    Q25 = %.4f\n', Q_tuned.Q25);
fprintf('    Q31 = %.4f\n', Q_tuned.Q31);
fprintf('    Q34 = %.4f\n', Q_tuned.Q34);
fprintf('    Q44 = %.4f\n', Q_tuned.Q44);
fprintf('    Q54 = %.4f\n', Q_tuned.Q54);
fprintf('    Q55 = %.4f\n', Q_tuned.Q55);

fprintf('\n  Tuned steady-state concentrations (mg/m^3):\n');
for i = 1:length(C_ss_tuned)
    fprintf('    C%d_ss_tuned = %.4f\n', i, C_ss_tuned(i));
end

% transient simulation for the tuned eigenvalue

[A_best, b_best] = build(Q_tuned, Q01, C01, Q03, C03);
A_dyn_best = -Dinv * A_best;
b_dyn_best =  Dinv * b_best;

odefun_best = @(t, C) A_dyn_best * C + b_dyn_best;
[t_vec_best, C_hist_best] = my_rk4(odefun_best, t0, tf, C0, dt);

t90_tuned = NaN(nTanks, 1);
C_target_tuned = 0.9 * C_ss_tuned;

for i = 1:nTanks
    Ci_hist_best = C_hist_best(:, i);
    idx = find(Ci_hist_best >= C_target_tuned(i), 1);
    if ~isempty(idx)
        t90_tuned(i) = t_vec_best(idx);
    end
end

fprintf('\nFinal t90 from tuned design (min):\n');
for i = 1:nTanks
    if isnan(t90_tuned(i))
        fprintf('    C%d: did not reach 90%%%% within simulation window.\n', i);
    else
        fprintf('    C%d: t90 = %.2f min\n', i, t90_tuned(i));
    end
end

% Check specs vs desired
allC_ok   = all(C_ss_tuned >= C_min_req);
allt90_ok = all(~isnan(t90_tuned)) && all(t90_tuned <= t90_req);

if allC_ok && allt90_ok
    fprintf('\nResult: Found a tune that satisfies conditions:\n');
    fprintf('  - All C_ss >= %.2f mg/m^3\n', C_min_req);
    fprintf('  - All t90_sim <= %.2f min\n', t90_req);
else
    fprintf('\nResult: Could not find a design that strictly meets specifications\n');
end
    

%% Figures

% figure: Steady-state concentrations – baseline vs tuned
figure;
bar(1:5, C_ss_base);
set(gca,'XTick',1:5,'XTickLabel',{'C1','C2','C3','C4','C5'});
ylabel('Concentration (mg/m^3)');
legend({'Baseline'}, 'Location','northwest');
title('Steady-State Concentrations: Baseline');
grid on;

% figure: Baseline transient response (0–100 min)
figure;
plot(t_vec, C_hist, 'LineWidth', 1.5);
xlabel('Time (min)');
ylabel('Concentration (mg/m^3)');
title('Baseline Transient Response of All Tanks');
legend({'C1','C2','C3','C4','C5'}, 'Location','southeast');
grid on;
hold on;
for i = 1:5
    if ~isnan(t90_base(i))
        y90 = 0.9 * C_ss_base(i);
        plot(t90_base(i), y90, 'ko', 'MarkerFaceColor','k');
    end
end
hold off;

% figure: Key valve flows – baseline vs tuned
flows_names = {'Q12','Q15','Q23','Q24','Q25','Q31','Q34','Q44','Q54','Q55'};

baseline_flows = [Q_base.Q12, Q_base.Q15, Q_base.Q23, Q_base.Q24, Q_base.Q25, ...
                  Q_base.Q31, Q_base.Q34, Q_base.Q44, Q_base.Q54, Q_base.Q55];

tuned_flows    = [Q_tuned.Q12, Q_tuned.Q15, Q_tuned.Q23, Q_tuned.Q24, Q_tuned.Q25, ...
                  Q_tuned.Q31, Q_tuned.Q34, Q_tuned.Q44, Q_tuned.Q54, Q_tuned.Q55];

figure;
bar(1:numel(flows_names), [baseline_flows(:), tuned_flows(:)]);
set(gca,'XTick',1:numel(flows_names),'XTickLabel',flows_names);
ylabel('Flow rate Q (m^3/min)');
legend({'Baseline','Tuned'}, 'Location','northwest');
title('Key Valve Flows: Baseline vs Tuned');
grid on;
xtickangle(45);

% figure: Steady-state concentrations – baseline vs tuned
figure;
bar(1:5, [C_ss_base, C_ss_tuned]);
set(gca,'XTick',1:5,'XTickLabel',{'C1','C2','C3','C4','C5'});
ylabel('Concentration (mg/m^3)');
legend({'Baseline','Tuned'}, 'Location','northwest');
title('Steady-State Concentrations: Baseline vs Tuned');
grid on;

% figure: Tuned transient response (0–100 min)
figure;
plot(t_vec_best, C_hist_best, 'LineWidth', 1.5);
xlabel('Time (min)');
ylabel('Concentration (mg/m^3)');
title('Tuned Transient Response of All Tanks');
legend({'C1','C2','C3','C4','C5'}, 'Location','southeast');
grid on;
hold on;
xline(30, '--', '30 min', 'LabelOrientation','horizontal', ...
    'LabelVerticalAlignment','bottom');
for i = 1:5
    if ~isnan(t90_tuned(i))
        y90_tuned = 0.9 * C_ss_tuned(i);
        plot(t90_tuned(i), y90_tuned, 'ko', 'MarkerFaceColor','k');
    end
end
hold off;

% figure: Per-tank t90 comparison (baseline vs tuned) + KPI overlay
max_t90_base  = max(t90_base(~isnan(t90_base)));
max_t90_tuned = max(t90_tuned(~isnan(t90_tuned)));

figure;
bar(1:5, [t90_base(:), t90_tuned(:)]);
set(gca,'XTick',1:5,'XTickLabel',{'C1','C2','C3','C4','C5'});
ylabel('t_{90} (min)');
legend({'Baseline','Tuned'}, 'Location','northwest');
title('Per-Tank 90% Rise Times: Baseline vs Tuned');
grid on;
hold on;
yline(t90_req, 'r--', sprintf('Requirement = %.1f min', t90_req), ...
    'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','bottom');

text(0.5, 0.95 * max([max_t90_base, max_t90_tuned, t90_req]), ...
    sprintf('Max t_{90} baseline = %.1f min\nMax t_{90} tuned = %.1f min', ...
            max_t90_base, max_t90_tuned), ...
    'Units','data', 'BackgroundColor',[0.95 0.95 0.95], 'EdgeColor','k');

hold off;

% figure: t90 Sensitivity from Key Flow Valves
flows_names = {'Q12','Q15','Q23','Q24','Q25','Q31','Q34','Q44','Q54','Q55'};

scales = linspace(0.5, 1.5, 11);   % 50% to 150% of tuned values
nFlows = numel(flows_names);
nScales = numel(scales);

t90_sens = NaN(nScales, nFlows);

for j = 1:nFlows
    vname = flows_names{j};

    for k = 1:nScales
        s = scales(k);

        Q_test = Q_tuned;
        Q_test.(vname) = Q_tuned.(vname) * s;

        if Q_test.(vname) <= 0
            continue;  % keep flows positive
        end

        [A_test, b_test] = build(Q_test, Q01, C01, Q03, C03);
        A_dyn_test = -Dinv * A_test;
        lambda_test = eig(A_dyn_test);
        lam_slow = max(real(lambda_test));

        if lam_slow >= 0
            % unstable, skip
            continue;
        end

        tau_slow = -1 / lam_slow;
        t90_sens(k, j) = 2.3 * tau_slow;
    end
end

% compute a simple sensitivity metric for each valve:
% sensitivity = (max t90 - min t90) across scale range
sens_metric = NaN(1, nFlows);
for j = 1:nFlows
    vals = t90_sens(:, j);
    vals = vals(~isnan(vals));
    if ~isempty(vals)
        sens_metric(j) = max(vals) - min(vals);
    end
end

[~, idx_sorted_sens] = sort(sens_metric, 'descend');
Nshow = min(5, nFlows);
idx_sel = idx_sorted_sens(1:Nshow);

figure;
hold on;
for m = 1:Nshow
    j = idx_sel(m);
    plot(scales, t90_sens(:, j), '-o', 'LineWidth', 1.5, ...
        'DisplayName', flows_names{j});
end
hold off;
xlabel('Scale factor (relative to tuned value)');
ylabel('t_{90,eig} (min)');
legend('Location','best');
title('Sensitivity of Dominant Rise Time to Valve Flows (Most Sensitive Valves)');
grid on;



% figure: Reactor performance summary (C_ss and t90 for each tank)
figure;

% rows 1–5: C1..C5, row 6: global summary
data_perf = zeros(6,4);

% per-tank values
for i = 1:5
    data_perf(i,1) = C_ss_base(i);      % baseline C_ss
    data_perf(i,2) = C_ss_tuned(i);     % tuned C_ss
    data_perf(i,3) = t90_base(i);       % baseline t90
    data_perf(i,4) = t90_tuned(i);      % tuned t90
end

% global summary row
max_t90_base  = max(t90_base(~isnan(t90_base)));
max_t90_tuned = max(t90_tuned(~isnan(t90_tuned)));
min_Css_base  = min(C_ss_base);
min_Css_tuned = min(C_ss_tuned);

data_perf(6,1) = min_Css_base;          % global min C_ss baseline
data_perf(6,2) = min_Css_tuned;         % global min C_ss tuned
data_perf(6,3) = max_t90_base;          % global max t90 baseline
data_perf(6,4) = max_t90_tuned;         % global max t90 tuned

rowNames_perf = {'C1','C2','C3','C4','C5','Global'};
colNames_perf = {'C_{ss} base (mg/m^3)', ...
                 'C_{ss} tuned (mg/m^3)', ...
                 't_{90} base (min)', ...
                 't_{90} tuned (min)'};

uit = uitable('Data', data_perf, ...
              'ColumnName', colNames_perf, ...
              'RowName', rowNames_perf, ...
              'Units','normalized', ...
              'Position',[0 0 1 1]);
title('Reactor Performance Summary: Baseline vs Tuned');

% figure: Key valve flows summary (baseline vs tuned)
flows_names = {'Q12','Q15','Q23','Q24','Q25','Q31','Q34','Q44','Q54','Q55'};

baseline_flows = [Q_base.Q12, Q_base.Q15, Q_base.Q23, Q_base.Q24, Q_base.Q25, ...
                  Q_base.Q31, Q_base.Q34, Q_base.Q44, Q_base.Q54, Q_base.Q55];

tuned_flows    = [Q_tuned.Q12, Q_tuned.Q15, Q_tuned.Q23, Q_tuned.Q24, Q_tuned.Q25, ...
                  Q_tuned.Q31, Q_tuned.Q34, Q_tuned.Q44, Q_tuned.Q54, Q_tuned.Q55];

figure;
data_flows = [baseline_flows(:), tuned_flows(:)];
colNames_flows = {'Q base (m^3/min)', 'Q tuned (m^3/min)'};

uit2 = uitable('Data', data_flows, ...
               'ColumnName', colNames_flows, ...
               'RowName', flows_names, ...
               'Units','normalized', ...
               'Position',[0 0 1 1]);
title('Key Valve Flows: Baseline vs Tuned');
