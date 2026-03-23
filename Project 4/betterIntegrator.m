function [I, counts] = betterIntegrator(x, y)
% Function: betterIntegrator.m
%
% This function numerically integrates a signal y over x using an adaptive
% approach combining trapezoidal, Simpson 1/3, and Simpson 3/8 rules.
% The method detects nonuniform intervals and applies trapezoidal rule
% where uniform spacing is not met. Counts of method usage are tracked.
%
% Inputs:
%   x      = independent variable vector (time or position)
%   y      = dependent variable vector (signal to integrate)
%
% Outputs:
%   I      = total integral of y over x
%   counts = vector of method usage counts: [trapezoidal, simpson 1/3, simpson 3/8]
%
% Author:       Isaac Chambers
% Section:      ME 2016 - C
% Assignment:   MP4
% Date:         11/23/2025

%% Initialize outputs and method counters
I = 0; 
counts = [0,0,0]';           % [trap, simpson 1/3, simpson 3/8] usage
trap_count = 0;
simp_one_third_count = 0;
simp_three_eigth_count = 0;

%% Loop through intervals, apply integration rules
i = 1;
total_intervals = length(x) - 1;

while i <= total_intervals - 3
    temp_x = x(i:i+2);       % check next 3 points for uniform spacing
    
    if ~isuniform(temp_x)    % nonuniform -> trapezoid
        I = I + trapz(temp_x(1:2), y(i:i+1));
        trap_count = trap_count + 1; 
        i = i + 1;
    else                     % uniform -> Simpson 1/3
        I_temp = (1/6)*(temp_x(3)-temp_x(1))*(y(i)+4*y(i+1)+y(i+2)); 
        I = I + I_temp; 
        simp_one_third_count = simp_one_third_count + 1; 
        i = i + 2;
    end
end

% Handle remaining intervals at end of dataset
remaining_intervals = total_intervals - (i-1);

if remaining_intervals == 0
    return;                 % no intervals left

elseif remaining_intervals == 1   % one interval left -> trapezoid
    I = I + trapz(x(end-1:end), y(end-1:end));
    trap_count = trap_count + 1;

elseif remaining_intervals == 2   % two intervals left
    temp_x = x(end-2:end);
    if isuniform(temp_x)
        I_temp = (1/6)*(temp_x(end)-temp_x(1))*(y(end-2)+4*y(end-1)+y(end));
        I = I + I_temp;
        simp_one_third_count = simp_one_third_count + 1;
    else
        I = I + trapz(x(end-2:end-1), y(end-2:end-1));
        I = I + trapz(x(end-1:end),   y(end-1:end));
        trap_count = trap_count + 2;
    end

elseif remaining_intervals == 3   % three intervals left
    temp_x = x(end-3:end);
    if isuniform(temp_x)    % uniform -> Simpson 3/8
        I_temp = (1/8)*(temp_x(end)-temp_x(1))*(y(end-3)+3*y(end-2)+3*y(end-1)+y(end));
        I = I + I_temp;
        simp_three_eigth_count = simp_three_eigth_count + 1;
    else                    % check subintervals for uniformity
        if isuniform(temp_x(1:3))
            I_temp = (1/6)*(x(end-1)-x(end-3))*(y(end-3)+4*y(end-2)+y(end-1));
            I = I + I_temp;
            simp_one_third_count = simp_one_third_count + 1;
            
            I = I + trapz(x(end-1:end), y(end-1:end));
            trap_count = trap_count + 1;
        else
            I = I + trapz(x(end-3:end-2), y(end-3:end-2));
            I = I + trapz(x(end-2:end-1), y(end-2:end-1));
            I = I + trapz(x(end-1:end),   y(end-1:end));
            trap_count = trap_count + 3;
        end
    end
end

% Return integral and method counts
counts = [trap_count, simp_one_third_count, simp_three_eigth_count];
