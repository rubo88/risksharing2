function ssr = griderror(coef)
global beta sigma b_grid nb  y_grid maxb minb eta

%% Parameters
    %nb =  size(b_grid,1);     % Size of capital grid
    %ny =  size(y_grid,1);     % Size of productivity grid
    y_L =y_grid(1);
    y_H =y_grid(2);
    N_p=int32(length(coef)/2);
    %q_number =  size(q_nodes,1);    % # of nodes
    p_c=coef(1:N_p);
    p_q=coef(N_p+1:end);
%% Sum erros over the grid   
c_L=NaN(nb,1);c_H=NaN(nb,1);
cnew_L=NaN(nb,1);cnew_H=NaN(nb,1);
q=NaN(nb,1);
b_prime=NaN(nb,1);
for i_b = 1:nb
    b           = b_grid(i_b);
    q(i_b)      =Phi(b/maxb,p_q);
    c_L(i_b)    =exp(Psi(b/maxb,p_c));
    c_L(i_b)    =min(min(c_L(i_b),y_L+b-q(i_b) *minb),y_H+y_L);
    c_H(i_b)    =y_L+y_H-c_L(i_b);
    c_H(i_b)    =min(c_H(i_b),y_H-b-q(i_b) *maxb);
    b_prime(i_b)=(y_L+b-c_L(i_b))./q(i_b);
    b_prime(i_b)=min(max(b_prime(i_b),minb),maxb);
    cnew_L(i_b)=exp(Psi(b_prime(i_b)/maxb,p_c));
    cnew_L(i_b)=min(min(cnew_L(i_b),y_L+b_prime(i_b)-q(i_b) *minb),y_H+y_L);
    cnew_H(i_b)=(y_L+y_H)-cnew_L(i_b);
    cnew_H(i_b)=min(cnew_H(i_b),y_H-b_prime(i_b)-q(i_b) *maxb);
end
%cnew_L=interp1(b_grid,c_L,b_prime,'linear','extrap');
%cnew_H=interp1(b_grid,c_H,b_prime,'linear','extrap');
% cnew_L=exp(Psi(b_prime/maxb,p_c))';
% cnew_H=(y_L+y_H).*ones(nb,1)-cnew_L;

ee_L   =q.*(c_L).^(-sigma)-beta.*eta.*(cnew_L).^(-sigma)-beta.*(1-eta).*(cnew_H).^(-sigma);
ee_H   =q.*(c_H).^(-sigma)-beta.*eta.*(cnew_H).^(-sigma)-beta.*(1-eta).*(cnew_L).^(-sigma);    
    
ssr=sum((ee_L).^2)+sum((ee_H).^2);

    
 

