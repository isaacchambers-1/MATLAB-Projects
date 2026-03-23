function [e_s] = significantDigitsPercentTolerance(n)
% This function calculates a conservative estimate for the relative percent
% specified tolerance (e_s) that gives accuracy to at least n significant
% digits. The formula is from Scarborough, 1996 / class notes. 
%
% Inputs:
%   n   =   integer, number of significant digits required
%
% Outputs:
%   e_s =   double, relative percent specified tolerance corresponding to n
%           significant digits
%
% Author:       Isaac Chambers, Graham Shepard
% Section:      ME 2016 - C
% Assignment:   MP1
% Date:         9/21/2025

    e_s = 0.5 * 10^(2-n);

end

