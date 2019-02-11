%
% very simple program to solve the standard growth model with projection
% methods
%
% simple grid (although ideally Chebyshev nodes should be used)
% regular polynomial (although ideally Cheb orth. pol should be used)
%
% the program uses two functions
% consfun(k,z,coef) evaluates consumption
% griderror(coef,k_grid,z_grid) calculates sum of squared Euler eq errors 

clear all
clc
global beta sigma b_grid nb Py ny y_grid maxb minb eta

%% parameter values
    ny = 2;               % number of states for the efficiency shock.
    beta   = 0.965;             % subjective discount factor
    y_grid=[0.9904 1.0470];
    y_L =y_grid(1);
    y_H =y_grid(2);
    Py  = [ .7412 , 0.2588 ;  0.2588, .7412];  
    eta=.7412 ;
    sigma = 3;                 % Coefficient of relative risk-aversion.
    minb =   -0.5;             % minimum value of the capital grid
    maxb =  0.5;                % maximum value of the capital grid 
    nb   = 101;               % number of grid points           
    b_grid = linspace(minb,maxb,nb)'; % grid equally spaced

%% initial value
    coef = [1 1];
    N_p=2;
%     load('initial')
%% find solution by minimizing the errors on the grid
for i=1:20
    options = optimset('MaxFunEvals',100000,'MaxIter',1000000,'TolFun',0.0000001,'TolX',0.0000001);
    coef_out = fminsearch(@mc,coef,options);
    coef=[coef_out 0];
    N_p=N_p+1;
    length(coef)
    mc(coef)
end

