function [y_fit, residuals] = lsRegression(x_data, y_data, fitType)

% Function: lsRegression.m
%
% This function performs least-squares regression for a given dataset and
% specified model type. It constructs the corresponding design matrix A 
% based on the regression form (linear, quadratic, cubic, power, logarithmic,
% or exponential), solves for the best-fit coefficients using the 
% pseudoinverse relation (AᵀA)⁻¹Aᵀb, and computes both the fitted model 
% predictions and residuals. 
%
% To handle datasets that include negative y-values, the function applies 
% a translation shift such that all y-values are nonnegative during fitting 
% (important for power and exponential models). The fitted results are then 
% shifted back to preserve the original scale.
%
% Supported model types:
%   'linear'           y = b₀ + b₁x
%   'quadratic'        y = b₀ + b₁x + b₂x²
%   'cubic'            y = b₀ + b₁x + b₂x² + b₃x³
%   'power'            y = a·xᵇ
%   'exponential'      y = a·e^(b·x)
%   'log'              y = b₀ + b₁·ln(x)
%
% Inputs:
%   x_data     =   vector of independent variable values (double)
%   y_data     =   vector of dependent variable values (double)
%   fitType    =   string specifying regression model type
%
% Outputs:
%   y_fit      =   fitted model output evaluated at x_data (double)
%   residuals  =   vector of y_data - y_fit (double)
%
% Author:       Isaac Chambers
% Section:      ME 2016 - C
% Assignment:   MP3
% Date:         10/29/2025


if (min(y_data) < 0)
    offset = abs(min(y_data(:)));
    y_translated = y_data(:) + offset;
else
    offset = 0;
    y_translated = y_data(:);
end

switch lower(fitType)

    case 'linear'
        A = [ones(size(x_data(:))), x_data(:)];
        b = pinv(A) * y_translated(:);
        y_fit = A * b;

    case 'quadratic'
        A = [ones(size(x_data(:))), x_data(:), x_data(:).^2];
        b = pinv(A) * y_translated(:);
        y_fit = A * b;

    case 'cubic'
        A = [ones(size(x_data(:))), x_data(:), x_data(:).^2, x_data(:).^3];
        b = pinv(A) * y_translated(:);
        y_fit = A * b;

    case 'power'
        A = [ones(size(x_data(:))), log(x_data(:) + eps)];
        b = pinv(A) * log(y_translated(:) + eps);
        y_fit = exp(b(1)) * (x_data(:) .^ b(2));

    case 'exponential'
        A = [ones(size(x_data(:))), x_data(:)];
        b = pinv(A)  * log(y_translated(:) + eps);
        y_fit = exp(b(1)) * exp(b(2) * x_data(:));

    case 'log'
        A = [ones(size(x_data(:))), log(x_data(:) + eps)];
        b = pinv(A) * y_translated(:);
        y_fit = A * b;

    otherwise
        error('Invalid fitType');

end

y_fit = y_fit - offset;
residuals = y_data(:) - y_fit;

end
