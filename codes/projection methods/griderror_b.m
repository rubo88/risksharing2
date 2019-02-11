function ssr = griderror_b(coef)
global beta sigma b_grid nb  y_grid maxb minb eta

%% Parameters
    %nb =  size(b_grid,1);     % Size of capital grid
    %ny =  size(y_grid,1);     % Size of productivity grid
    yL =y_grid(1);
    yH =y_grid(2);
    N_p=int32(length(coef)/2);
    %q_number =  size(q_nodes,1);    % # of nodes
    p_b=coef(1:N_p);
    p_q=coef(N_p+1:end);
%% Sum erros over the grid   
c_L=NaN(nb,1);c_H=NaN(nb,1);
c_LL=NaN(nb,1);c_LH=NaN(nb,1);c_HL=NaN(nb,1);c_HH=NaN(nb,1);
q=NaN(nb,1);
bp_L=NaN(nb,1);bp_H=NaN(nb,1);
bp_LL=NaN(nb,1);bp_HH=NaN(nb,1);
bp_LH=NaN(nb,1);bp_HL=NaN(nb,1);
for i_b = 1:nb
    b           = b_grid(i_b);
    q(i_b)      =exp(Phi(b/maxb,p_q));
    %TODAY
        bp_L(i_b) =Psi2(b/maxb,p_b)*maxb;    
        c_L(i_b)  =yL+b-q(i_b)*bp_L(i_b);
             if c_L(i_b)<=0
                c_L(i_b)    =0;
                bp_L(i_b) =(yL+b)/q(i_b);
            elseif c_L(i_b)>=yH+yL
                c_L(i_b)=yH+yL;
                bp_L(i_b) =(-yH+b)/q(i_b);
             end
        bp_H(i_b) =-bp_L(i_b);
        c_H(i_b)  =yH+yL-c_L(i_b);
    %TOMORROW
        q2          =exp(Phi(bp_L(i_b)/maxb,p_q));
    %LL
        bp_LL(i_b)=Psi2(bp_L(i_b)/maxb,p_b)*maxb;
        c_LL(i_b) =yL+bp_L(i_b)-q2*bp_LL(i_b);
            if c_LL(i_b)<=0
                c_LL(i_b)    =0;
                bp_LL(i_b) =(yL+bp_L(i_b))/q(i_b);
            elseif c_LL(i_b)>=yH+yL
                c_LL(i_b)=yH+yL;
                bp_LL(i_b) =(-yH+bp_L(i_b))/q(i_b);
            end
        bp_HH(i_b)=-bp_LL(i_b);    
        c_HH(i_b)=(yL+yH)-c_LL(i_b);    
    %LH
        bp_HL(i_b)=Psi2(bp_H(i_b)/maxb,p_b)*maxb;
        c_HL(i_b)   =yL+bp_H(i_b)-q2*bp_HL(i_b);
            if c_HL(i_b)<=0
                c_HL(i_b)    =0;
                bp_HL(i_b) =(yL+bp_H(i_b))/q2;
            elseif c_HL(i_b)>=yH+yL
                c_HL(i_b)=yH+yL;
                bp_HL(i_b) =(-yH+bp_H(i_b))/q2;
            end 
        bp_LH(i_b)=-bp_HL(i_b);
        c_LH(i_b)=(yL+yH)-c_HL(i_b);
end

ee_L   =q.*(c_L).^(-sigma)-beta.*eta.*(c_LL).^(-sigma)-beta.*(1-eta).*(c_LH).^(-sigma);
ee_H   =q.*(c_H).^(-sigma)-beta.*eta.*(c_HH).^(-sigma)-beta.*(1-eta).*(c_HL).^(-sigma);    
    
ssr=sum((ee_L).^2)+sum((ee_H).^2);

    
 

