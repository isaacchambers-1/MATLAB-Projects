function [ x_r, iter, time ] = newtonRaphsonMethod( f, dfdx, x_0, e_s )

% This function applies the Newton–Raphson iterative root-finding method to
% locate a zero of the function f(x), using its derivative dfdx(x). The method
% begins from an initial guess and updates the estimate using tangent-line
% intersections until the relative error falls below the specified tolerance or a
% maximum iteration count is reached. The function returns the final root
% estimate, total iterations used, and elapsed execution time.

% Inputs:
%   f       =   function handle representing f(x)
%   dfdx    =   function handle representing the derivative f'(x) 
%   x_0     =   initial guess for the root (double)
%   e_s     =   stopping tolerance based on approximate relative error (%)

% Outputs:
%   x_r     =   final root estimate obtained from Newton–Raphson
%   iter    =   total number of iterations executed before convergence
%   time    =   total runtime for the method (seconds), measured via tic/toc

% Author:       Isaac Chambers, Graham Shepard
% Section:      ME 2016 - C
% Assignment:   MP2
% Date:         10/5/2025

tic % begin performance timer

rel_error = 100; % initialize high error to enter the loop
iteration = 0; 
max_iteration = 1000; % prevent the loop from repeating infinitely in the case of divergence 

% inital check if guess is a root itself
if f(x_0) == 0
    x_r = x_0;
    iter = iteration; 
    time = toc; 
    return; 
end 

% derivates very close to 0 cause instability and can sometimes cause
% divergence - warns user of such 
while rel_error > e_s && iteration < max_iteration
    if abs(dfdx(x_0)) < 1e-12
        warning("Derivative is near zero. Method may fail.")
        break; 
    end
    
    x_next = x_0 - (f(x_0) / dfdx(x_0)); % core newton raphson formula implemented
    rel_error = abs((x_next - x_0) / x_next) * 100; % use relative approximate error to determine convergence according to inputted e_s
    x_0 = x_next; 
    iteration = iteration + 1; 

end

if iteration == max_iteration
    error('Newton-Raphson did not converge within the maximum iteration limit.');
end

x_r = x_0; 
iter = iteration; 
time = toc;

