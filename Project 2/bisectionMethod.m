function [ x_r, iter, time ] = bisectionMethod( f, x_l, x_u, e_s )

% This function applies the Bisection Method to locate a root of the function f(x)
% within a specified interval [x_l, x_u]. It assumes the function changes sign
% across the interval (f(x_l)*f(x_u) < 0), bracketing a root, and repeatedly
% halves the interval on the proper side where there is still a sign change, 
% and therefore still surrounds the root, until the approximate relative error
% is acceptable. The function returns the estimated root, number of iterations, 
% and total execution time.

% Inputs:
%   f       =   function handle representing f(x)
%   x_l     =   lower bound of bracketing interval (double)
%   x_u     =   upper bound of bracketing interval (double)
%   e_s     =   stopping tolerance based on relative approximate error (%)

% Outputs:
%   x_r     =   final root estimate obtained from bisection
%   iter    =   total number of iterations executed
%   time    =   elapsed runtime in seconds, measured using tic/toc

% Author:       Isaac Chambers, Graham Shepard
% Section:      ME 2016 - C
% Assignment:   MP2
% Date:         10/5/2025

tic % begin performance timer

relative_error = 100; 
iteration = 0; 

% validate that the inital guess actually brackets a root 
if (f(x_l) * f(x_u) > 0)
    error('Initial guess invalid. See help doc for bisectionMethod.')
end

% check that the root doesn't fall exactly on an endpoint
if f(x_l) == 0
    x_r = x_l; 
    iter = 0; 
    time = toc; 
    return;
    
elseif f(x_u) == 0
    x_r = x_u; 
    iter = 0; 
    time = toc; 
    return;
end


% bisection method - continue halving the interval on the right side so
% the new interval brackets the root until the desired tolerance is met
while (relative_error > e_s)
    x_root = .5 * (x_l+x_u); % midpoint as current guess
    if (f(x_root) == 0 ) % just in case the exact root is landed on
        x_r = x_root;
        iter = iteration;
        time = toc;
    elseif (f(x_l) * f(x_root) < 0) % root is in left half 
        x_l = x_l; %redudant update for clarity 
        x_u = x_root;
    elseif (f(x_l) * f(x_root) > 0) % root is in right half 
        x_l = x_root;
        x_u = x_u; % redundant update for clarity
    else
        error('Something went wrong')
    end 

    x_r_next = .5 * (x_l+x_u); 
    relative_error = abs((x_r_next-x_root)/x_r_next) * 100; 
    iteration = iteration + 1; 
end

x_r = x_r_next;
iter = iteration;
time = toc;

end


