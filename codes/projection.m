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
global beta sigma b_grid nb y_grid maxb minb eta

%% parameter values
    ny = 2;               % number of states for the efficiency shock.
    beta   = 0.965;             % subjective discount factor
    y_grid=[0.9904 1.0470];
    y_L =y_grid(1);
    y_H =y_grid(2);
    Py  = [ .7412 , 0.2588 ;  0.2588, .7412];  
    eta=.7412 ;
    sigma = 3;                 % Coefficient of relative risk-aversion.
    minb =   -1;             % minimum value of the capital grid
    maxb =  1;                % maximum value of the capital grid 
    nb   = 101;               % number of grid points           
    b_grid = linspace(minb,maxb,nb)'; % grid equally spaced

%% initial value
    coef_c = [1 -1];
    N_p=length(coef_c);
    coef_q = [0 1];
    coef_in= [coef_c coef_q];
    
%     coef_c = [1.3574   -1.5384   -0.0027    0.0001    0.2175    0.0026   -0.0199    0.0017    0.0018   -0.0012   -0.0013         0];
%     N_p=length(coef_c);
%     coef_q = [0.0007   -0.0013    0.0003   -0.0000   -0.0003    0.0000   -0.0000   -0.0001   -0.0000   -0.0000   -0.0000         0];
%     coef_in= [coef_c coef_q];
    
    %load('initial')
%% find solution by minimizing the errors on the grid
for i=1:20
    options = optimset('MaxFunEvals',100000,'MaxIter',1000000,'TolFun',0.0000001,'TolX',0.0000001);
    coef_out = fminsearch(@griderror,coef_in,options);
    coef_c_out=[coef_out(1:N_p) 0];
    coef_q_out=[coef_out(N_p+1:end) 0];
    N_p=N_p+1;
    coef_in=[coef_c_out coef_q_out];
    i
    griderror(coef_out)
end
%% plot the consumption choice as a function of k (for 3 values of z)
figure
c_L=exp(Psi(b_grid./maxb,coef_c_out)-1)*(y_H+y_L);
plot(b_grid,c_L);
figure
q=Phi(b_grid./maxb,coef_q_out);
plot(b_grid,q);
figure
plot(b_grid,(y_L+b_grid-c_L')./q');
