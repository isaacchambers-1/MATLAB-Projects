function [t_vec, C_hist] = my_rk4(odefun, t0, tf, C0, dt)
% my_rk4
%
% This function implements a fixed-step 4th-order Runge-Kutta (RK4) method
% to numerically integrate a system of first-order ODEs of the form
%   dC/dt = odefun(t, C)
% over the time interval [t0, tf]. It returns the time vector and the
% history of the state vector at each time step.
%
% Inputs:
%   odefun  = function handle for the ODE right-hand side:
%             dCdt = odefun(t, C)
%   t0      = initial time (min)
%   tf      = final time (min)
%   C0      = initial state vector at t0 (can be row or column)
%   dt      = fixed time step size (min)
%
% Outputs:
%   t_vec   = column vector of times [t0, t0+dt, ..., t0+N*dt] (min)
%   C_hist  = matrix of state values; each row corresponds to the state
%             vector at the corresponding time in t_vec
%
% Author: Isaac CHambers
% Section: ME 2016
% Assignment: MP5
% Date: 12/11/25

  

    % Number of time steps (integer count based on dt)
    n_steps = floor((tf - t0) / dt);

    % Time vector (column): uniform grid from t0 in increments of dt
    t_vec = (0:n_steps)' * dt + t0;

    % Ensure initial condition is a column vector
    C0 = C0(:);
    n_states = length(C0);

    % Preallocate history matrix for efficiency
    C_hist = zeros(n_steps + 1, n_states);

    % Store initial condition as first row of history
    C_hist(1, :) = C0.';   % transpose to row

   
    for k = 1:n_steps
        t_k = t_vec(k);
        C_k = C_hist(k, :).';   % current state as column vector

        % RK4 stages: sample slope at four locations within the step
        k1 = odefun(t_k,            C_k);
        k2 = odefun(t_k + dt/2.0,   C_k + dt/2.0 * k1);
        k3 = odefun(t_k + dt/2.0,   C_k + dt/2.0 * k2);
        k4 = odefun(t_k + dt,       C_k + dt      * k3);

        % Combine stages to obtain next state
        C_next = C_k + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4);

        % Store updated state in history (as row)
        C_hist(k + 1, :) = C_next.';
    end
end
