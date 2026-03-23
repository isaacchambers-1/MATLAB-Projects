function [D, counts] = betterDifferentiator(x, y)
% Function: betterDifferentiator.m
%
% Numerically differentiates y with respect to x using an adaptive approach
% that combines central difference (CDD) for uniform intervals and
% interpolating polynomial divided difference (IPDD) for nonuniform intervals.
% Counts of method usage are tracked.
%
% Inputs:
%   x      = independent variable vector (time or position)
%   y      = dependent variable vector (signal to differentiate)
%
% Outputs:
%   D      = derivative of y with respect to x
%   counts = vector of method usage counts: [CDD, IPDD]
%
% Author:       Isaac Chambers
% Section:      ME 2016 - C
% Assignment:   MP4
% Date:         11/23/2025

% Initialize outputs and counters
D = zeros(size(x));
counts = [0, 0]';          % [CDD, IPDD]
cdd_count = 0;
ipdd_count = 0;

% First point - forward difference or interpolating difference
if ~isuniform(x(1:3))        % nonuniform -> IPDD
    D(1) = y(1)*(2*x(1)-x(2)-x(3))/((x(1)-x(2))*(x(1)-x(3))) + ...
           y(2)*(2*x(1)-x(1)-x(3))/((x(2)-x(1))*(x(2)-x(3))) + ...
           y(3)*(2*x(1)-x(1)-x(2))/((x(3)-x(1))*(x(3)-x(2)));
    ipdd_count = ipdd_count + 1;
else                         % uniform -> forward CDD
    D(1) = (-3*y(1) + 4*y(2) - y(3)) / (2*(x(2)-x(1)));
    cdd_count = cdd_count + 1;
end

% Interior points - central or interpolating difference
i = 2;
while i <= length(x)-1
    temp_x = x(i-1:i+1);      % take current, previous, next points
    
    if ~isuniform(temp_x)     % nonuniform -> IPDD
        D_temp = y(i-1)*(2*x(i)-x(i)-x(i+1))/((x(i-1)-x(i))*(x(i-1)-x(i+1))) + ...
                 y(i)*(2*x(i)-x(i-1)-x(i+1))/((x(i)-x(i-1))*(x(i)-x(i+1))) + ...
                 y(i+1)*(2*x(i)-x(i-1)-x(i))/((x(i+1)-x(i-1))*(x(i+1)-x(i)));
        D(i) = D_temp;
        ipdd_count = ipdd_count + 1;
    else                      % uniform -> central difference
        D_temp = (y(i+1)-y(i-1)) / (2*(x(i+1)-x(i)));
        D(i) = D_temp;
        cdd_count = cdd_count + 1;
    end
    i = i + 1;
end

% Last point - backward difference or interpolating difference
if ~isuniform(x(end-2:end))   % nonuniform -> IPDD
    D(end) = y(end-2)*(2*x(end)-x(end-1)-x(end))/((x(end-2)-x(end-1))*(x(end-2)-x(end))) + ...
             y(end-1)*(2*x(end)-x(end-2)-x(end))/((x(end-1)-x(end-2))*(x(end-1)-x(end))) + ...
             y(end)*(2*x(end)-x(end-2)-x(end-1))/((x(end)-x(end-2))*(x(end)-x(end-1)));
    ipdd_count = ipdd_count + 1;
else                          % uniform -> backward CDD
    D(end) = (3*y(end) - 4*y(end-1) + y(end-2)) / (2*(x(end)-x(end-1)));
    cdd_count = cdd_count + 1;
end

% Return method counts
counts = [cdd_count, ipdd_count];
