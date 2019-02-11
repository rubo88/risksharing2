
clear all
clc
%global beta sigma b_grid nb y_grid maxb minb eta

%% parameter values
    ny = 2;               % number of states for the efficiency shock.
    beta   = 0.965;             % subjective discount factor
    y=[0.9904 1.0470];
    yL =y(1);
    yH =y(2);
    Py  = [ .7412 , 0.2588 ;  0.2588, .7412];  
    eta=.7412 ;
    gamma = 3;                 % Coefficient of relative risk-aversion.
    minb =   -1;             % minimum value of the capital grid
    maxb =  1;                % maximum value of the capital grid 
    nb   = 101;               % number of grid points           
    bgrid = linspace(minb,maxb,nb)'; % grid equally spaced
    
% My guess of cL policy:
    cL=yL+bgrid;cL(cL<=0)=0.0001;
    cH=yH-bgrid;cL(cH<=0)=0.0001;
    
stop_c=0;
%while stop_c==0
    % Guess a price function
    q=ones(nb,1);
    
    
    
    
    stop_q=0;i=0;
    while stop_q==0
        i=i+1
    % RHS of EE
    rhsL=beta*eta*(cL).^(-gamma)+beta*(1-eta)*(cH).^(-gamma);
    rhsH=beta*eta*(cH).^(-gamma)+beta*(1-eta)*(cL).^(-gamma);
    % C today on tomorrows grid
    cL_today=(rhsL./q).^(-1/gamma);
    cH_today=(rhsH./q).^(-1/gamma);
    % Implied b today
    bL_today=cL_today+q.*bgrid-yL;
    bH_today=cH_today+q.*bgrid-yH;
    % Interpolate c on the original grid
    cL_new=interp1(bL_today,cL_today,bgrid,'linear','extrap');
    cH_new=interp1(bH_today,cH_today,bgrid,'linear','extrap');
    % Check borrowing constraint
    bindingL=bgrid<repmat(bL_today(1),nb,1);bindingH=bgrid<repmat(bH_today(1),nb,1);
    AUXL=yL+bgrid-q.*minb;AUXH=yH+bgrid-q.*minb;
    cL_new(bindingL) = AUXL(bindingL);cH_new(bindingH) = AUXH(bindingH);
    
    qnew=q-(cL_new+cH_new>yH+yL).*(0.001)+(cL_new+cH_new<yH+yL).*(0.001);
    %qnew=q-(bL_today>bH_today).*(0.01)+(bL_today<bH_today).*(0.01);
    if sum(abs(qnew-q))<0.1
        stop_q=1;
    end
    plot(q)
    pause(0.1)
    q=qnew;
    end
    
    cL=cLnew;cH=cHnew;
    
    
    
%end
    
    
    
    
    