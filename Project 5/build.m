function [A, b] = build(Q, Q01, C01, Q03, C03)
% build
%
% Constructs the steady-state linear system A*C_ss = b for a given set of
% flow rates Q.* and source terms into Tanks 1 and 3. This is SPECIFICALLY
% for the Webb Industries current reactor tank setup. 
%
% Inputs:
%   Q         = structure containing all inter-tank flow rates
%   Q01, C01  = source flow and concentration into Tank 1
%   Q03, C03  = source flow and concentration into Tank 3
%
% Outputs:
%   A = 5x5 steady-state coefficient matrix
%   b = right-hand-side vector from inflow sources
%
% Author: Isaac Chambers
% Section: ME 2016 - C
% Assignment: MP5
% Date: 12/11/25

A = [ (Q.Q12 + Q.Q15),          0,                -Q.Q31,             0,               0;
      -Q.Q12,                   (Q.Q23 + Q.Q24 + Q.Q25),   0,        0,               0;
       0,                       -Q.Q23,           (Q.Q31 + Q.Q34),     0,              0;
       0,                       -Q.Q24,           -Q.Q34,             Q.Q44,          -Q.Q54;
      -Q.Q15,                   -Q.Q25,            0,                 0,     (Q.Q55 + Q.Q54)];

b = [Q01 * C01;
     0;
     Q03 * C03;
     0;
     0];

end
