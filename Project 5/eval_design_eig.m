function [cost_val, C_ss_out, t90_eig, lambda_slowest] = eval_design_eig(Q, Q01, C01, Q03, C03, Dinv, C_min_req, Q_nom)
% eval_design_eig
%
% Evaluates a proposed reactor-network tune by computing the steady-state
% concentrations, the dominant eigenvalue of the linearized dynamic system,
% and an associated cost based on settling time and constraint violations.
%
% Inputs:
%   Q         = structure containing design flow rates
%   Q01, C01  = source flow and concentration into Tank 1
%   Q03, C03  = source flow and concentration into Tank 3
%   Dinv      = diagonal matrix of 1/Vi terms
%   C_min_req = minimum allowable concentration across tanks
%   Q_nom     = nominal flow-rate structure 
%
% Outputs:
%   cost_val       = objective function value for optimization
%   C_ss_out       = steady-state concentration vector
%   t90_eig        = approximate t90 computed from slowest eigenvalue
%   lambda_slowest = eigenvalue governing slow system dynamics
%
% Author: Isaac Chambers
% Section: ME 2016 - C
% Assignment: MP5
% Date: 12/11/25

    wC = 20;
    wQ = 1.0;

    % Build linear steady-state system
    [A, b] = build(Q, Q01, C01, Q03, C03);

    % Reject numerically unusable designs
    if rcond(A) < 1e-10
        cost_val       = Inf;
        C_ss_out       = nan(5,1);
        t90_eig        = Inf;
        lambda_slowest = Inf;
        return
    end

    % Solve for steady-state concentrations
    C_ss_out = A \ b;

    % Enforce physical feasibility
    if any(C_ss_out < 0)
        cost_val       = Inf;
        t90_eig        = Inf;
        lambda_slowest = Inf;
        return
    end

    % Eigenvalue analysis
    A_dyn = -Dinv * A;
    lambda = eig(A_dyn);
    lambda_slowest = max(real(lambda));

    % Reject unstable or non-decaying designs
    if lambda_slowest >= 0
        cost_val       = Inf;
        t90_eig        = Inf;
        return
    end

    % Convert slowest eigenvalue to approximate t90
    tau     = -1 / lambda_slowest;
    t90_eig = 2.3 * tau;

    % Concentration constraint violation
    minC      = min(C_ss_out);
    violation = max(0, C_min_req - minC);

    Q_vec      = [Q.Q12; Q.Q15; Q.Q23; Q.Q24; Q.Q25; Q.Q31; Q.Q34; Q.Q44; Q.Q54; Q.Q55];
    Q_nom_vec  = [Q_nom.Q12; Q_nom.Q15; Q_nom.Q23; Q_nom.Q24; Q_nom.Q25; Q_nom.Q31; Q_nom.Q34; Q_nom.Q44; Q_nom.Q54; Q_nom.Q55];
    ratio_diff_sq = mean((Q_vec ./ Q_nom_vec - 1).^2);

    % Total cost
    cost_val = t90_eig + wC*violation + wQ*ratio_diff_sq;

end
