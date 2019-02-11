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
c_LL=NaN(nb,1);c_LH=NaN(nb,1);c_HL=NaN(nb,1);c_HH=NaN(nb,1);
q=NaN(nb,1);
bp_L=NaN(nb,1);bp_H=NaN(nb,1);
for i_b = 1:nb
    b           = b_grid(i_b);
    q(i_b)      =exp(Phi(b/maxb,p_q));
    
    c_L(i_b)    =exp(Psi(b/maxb,p_c)-1)*(y_H+y_L);
    c_H(i_b)    =y_L+y_H-c_L(i_b);
    
    bp_L(i_b)=(y_L+b-c_L(i_b))./q(i_b);
    bp_L(i_b)=(bp_L(i_b)<minb)*minb+(bp_L(i_b)>maxb)*maxb+(bp_L(i_b)<minb)*(bp_L(i_b)>maxb)*bp_L(i_b);
    
    c_LL(i_b)=exp(Psi(bp_L(i_b)/maxb,p_c)-1)*(y_H+y_L);
    c_HH(i_b)=(y_L+y_H)-c_LL(i_b);
    
    bp_H(i_b)=-bp_L(i_b);
    
    c_HL(i_b)=exp(Psi(bp_H(i_b)/maxb,p_c)-1)*(y_H+y_L);
    c_LH(i_b)=(y_L+y_H)-c_HL(i_b);
end


ee_L   =q.*(c_L).^(-sigma)-beta.*eta.*(c_LL).^(-sigma)-beta.*(1-eta).*(c_LH).^(-sigma);
ee_H   =q.*(c_H).^(-sigma)-beta.*eta.*(c_HH).^(-sigma)-beta.*(1-eta).*(c_HL).^(-sigma);    
    
ssr=sum((ee_L).^2)+sum((ee_H).^2);

    
 

