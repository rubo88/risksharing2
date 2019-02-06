% =====================================================================
% This program solves the Aiyagari self-insurance model ===============
% =====================================================================
clear;
clf;
tic
disp('Solving the Aiyagari self-insurance model by value function iterations')
disp('    ')
disp('Start by finding the equilibrium interest rate with a bisection method')

% ======================================================================
% ================        User Definition Area       ===================
% ======================================================================

nstates = 2;               % number of states for the efficiency shock.
beta   = 0.965;             % subjective discount factor
y=[0.9904 1.0470];
prob   = [ .7412 , 0.2588 ; .7412 , 0.2588];  
                           % prob(i,j) = probability (A(t+1)=Aj | A(t) =
                           % Ai) - THIS NEEDS TO BE CONSISTENT WITH NSTATES

sigma = 3;                 % Coefficient of relative risk-aversion.

minb =   -25;             % minimum value of the capital grid
maxb =  25;                % maximum value of the capital grid 
nb   = 1001;               % number of grid points             
db=0.001;
dq = 1;                    % sufficiently big initial difference of old and new r
rcrit = 1e-3;              % Convergence criteria for r
iterout = 0;               % initialize counter of iterations to 0     
tolVFI=1e-3;               % Convergence criteria for VFI

% ======================================================================
% =========   Ergodic dist, set initial r, create grids etc  ===========       
% ======================================================================
% Find the ergodic probabilities of the transition probability matrix
plr = prob^500;

% Following Aiyagari (1994), the initial interest rate is set as 1/beta-1-epsilon.
epsilon = 1e-5;
rguess=ones(nb,1).*(1/beta-1-epsilon);
qguess=1./(rguess);
% the grid is now computed               
bgrid = linspace(minb,maxb,nb)'; % grid equally spaced

% ======================================================================
% =========   THE MAIN LOOP STARTS HERE  ===============================       
% ======================================================================

q=ones(nb,1).*0.93;
q_new=zeros(nb,1);
while dq>rcrit
    iterout = iterout+1;    % Count the number of iterations

    % Form single period return function  
    % This is a hyper-matrix of dimension(nstates,nk,nk)
    % 1st dim is k', 2nd is k and 3rd is shock
    % In principle we set all values at -inf (you will see later why)
        RR=-10000000*ones(nb,nb,nstates); 
        
    %  loop through all possible states (i.e. K) 
    consump = zeros(nstates,1);
    for i=1:nb
        b=bgrid(i);
            %  loop through all possible controls (i.e. K')    
            for j=1:nb
            bprime = bgrid(j);
                % loop through all shocks
                for l=1:nstates
                    % Compute consump for every possible state
                    consump = y(l)+b-q(i)*bprime;
                    % Fill up matrix with utility for every state/control
                    if consump > 0
                        if sigma==1
                            RR(j,i,l) = log(consump);
                        else
                            RR(j,i,l) = (consump).^(1-sigma)./(1-sigma);
                        end
                    end
                    % Notice that if consump <0, utility is defined as
                    % a very negative number
                end
            end
     end

    % initialize some variables
        v       = zeros(nb,nstates);    % values
        tv      = zeros(nb,nstates);    % new iteration values
        decis   = zeros(nb,nstates);    % optimal policy for k' 
        tdecis  = zeros(nb,nstates);    % new optimal policy for k'
            % decis and tdecis tell you which is the optimal k' for every (s,k)
            % It doesn't tell you the value of k', but the position in the
            % grid! (therefore its filled with natural numbers)
        iter    = 0;
        metric1 = 1;metric2 = 1;                        
        metric = max(metric1,metric2);  % distances between two iterations     
      
% ======================================================================
% =========   VALUE FUNCTION ITERATION  ================================       
% ======================================================================
    %  Iterate on Bellman's equation and get the decision 
    %  rules and the value function at the optimum

        while metric>tolVFI
            % Get the optimal policy for every posible state
            for l=1:nstates                             % Loop on the shock
                RRR=zeros(nb,nb);               
                RRR(:,:)=RR(:,:,l);                     % Take the U(k,k') matrix
                zz=RRR+beta*(v*prob(l,:)')*ones(1,nb);  % This is Bellman equation
                [tv(:,l),tdecis(:,l)] = max(zz);        % Maximize for every k.
            end            
            % Check convergence (both in policy and value function)
                metric1 = max(abs(tdecis-decis));
                metric2 = max(abs(v-tv)./(1+abs(v)));
                metric = max(metric1,metric2);
            % Update value and policy functions
                v=tv;
                decis=tdecis;
            iter=iter+1;
        end
        
    % Compute the other optimal decisitions from the decision of k'
        bdecis = zeros(nb,nstates); 
        cdecis = zeros(nb,nstates);
        %idecis = zeros(nb,nstates);
        %indecis = zeros(nb,nstates);
        for i1=1:size(tdecis,1)
            for i2=1:size(tdecis,2)
                bdecis(i1,i2)  = bgrid(tdecis(i1,i2));
                    % Notice that kdecis containts the same information as
                    % tdecis (see above comment), but now it tells you the real
                    % value of i2' (a real number) that is optimal, 
                    % not the position of the grid
                cdecis(i1,i2)  = y(i2)+bgrid(i1)-q(i1)*bdecis(i1,i2);
                %idecis(i1,i2)  = bdecis(i1,i2)- (1-delta)*bgrid(i1);
                %indecis(i1,i2) = idecis(i1,i2) - delta*bgrid(i1);
            end
        end


% ======================================================================
% =========   NEW INTEREST RATE  =======================================
% ======================================================================
    for i=1:(nb-1)/2+1
        ed=bdecis(nb+1-i,2)+bdecis(i,1)
        if ed>db
            q_new(nb+1-i)=q(nb+1-i)+0.0001*abs(ed);
            q_new(i)=q_new(nb+1-i);
        elseif ed<-db
            q_new(nb+1-i)=q(nb+1-i)-0.0001*abs(ed);
            q_new(i)=q_new(nb+1-i);
        else
            q_new(nb+1-i)=q(nb+1-i);
            q_new(i)=q_new(nb+1-i);
        end
    end
    hold all
    plot(bgrid,q_new);
    pause(0.1)
    dq=abs(sum(q_new-q));
    if dq>rcrit
%         disp(['Initial r' '   ' 'Implied r'])
%         disp([r rnew])
        q = (q_new+q)/2;
    end

end
% ======================================================================
% =========   MAIN LOOP ENDS HERE  =====================================       
% ======================================================================

disp('Market clearing interest rate is found at:')
disp(r)
disp('   ')


%   simulate life histories of the agent

%   print out results 

disp('Value Function Iterations and Computation of Stationary Distribution and Simulation are Completed.  Time Spent:')
toc
disp('   ')
disp('PARAMETER VALUES')
disp('   ')
disp('    alfa      beta      delta  ')
disp([alfa, beta, delta])
disp('   ')
disp('Values of Productivity')
disp(aval)
disp('Transition Probability Matrix for A')
disp(prob)
disp('   ')
disp('RESULTS')
disp('   ')
disp('    mean of K ')
disp(         meanK)
disp(' as compared to mean of deterministic steady-states: ')
disp(steadyk)


plot(bgrid,bdecis);
title('Incomplete Markets Model: POLICY FUNCTION for kprime');
hold on
plot(bgrid,bgrid,'k');
hold off

pause

plot(bgrid,cdecis);
title('Incomplete Markets Model: POLICY FUNCTION for C');

pause

plot(bgrid,idecis);
title('Incomplete Markets Model: POLICY FUNCTION for I');

pause

plot(bgrid,probk');
title('Stationary DISTRIBUTION OF CAPITAL');

pause






