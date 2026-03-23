function [ rmseVal ] = rmsError( y_data, y_fit )

% Function: rmsError.m
%
% This function computes the Root Mean Square Error (RMSE) between a set of 
% observed (actual) data points and corresponding model predictions. RMSE 
% provides a quantitative measure of the model's overall predictive accuracy 
% by representing the square root of the mean of squared residuals.
%
% A lower RMSE value indicates a model that more closely fits the data. This 
% metric is especially useful when comparing the performance of multiple 
% regression models on the same dataset.
%
% Inputs:
%   y_data     =   vector of observed or measured data values (double)
%   y_fit      =   vector of predicted or fitted model values (double)
%
% Outputs:
%   rmseVal    =   scalar value of root mean square error (double)
%
% Author:       Isaac Chambers
% Section:      ME 2016 - C
% Assignment:   MP3
% Date:         10/29/2025


rmseVal = 0;

for i = 1:numel(y_data)
    rmseVal = rmseVal + (y_data(i) - y_fit(i))^2;
end

rmseVal = sqrt(rmseVal / numel(y_data));


