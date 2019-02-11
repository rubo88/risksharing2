function ssr=mc(coef)
global beta sigma b_grid nb Py ny y_grid maxb minb eta

 
tolVFI=1e-3;               % Convergence criteria for VFI

q=exp(Phi(b_grid./maxb,coef));

    % Form single period return function  
    % This is a hyper-matrix of dimension(nstates,nk,nk)
    % 1st dim is k', 2nd is k and 3rd is shock
    % In principle we set all values at -inf (you will see later why)
        RR=-10000000*ones(nb,nb,ny); 
        
    %  loop through all possible states (i.e. K) 
    consump = zeros(ny,1);
    for i=1:nb
        b=b_grid(i);
            %  loop through all possible controls (i.e. K')    
            for j=1:nb
            bprime = b_grid(j);
                % loop through all shocks
                for l=1:ny
                    % Compute consump for every possible state
                        consump = y_grid(l)+b-q(i)*bprime;
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
        v       = zeros(nb,ny);    % values
        tv      = zeros(nb,ny);    % new iteration values
        decis   = zeros(nb,ny);    % optimal policy for k' 
        tdecis  = zeros(nb,ny);    % new optimal policy for k'
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
            for l=1:ny                             % Loop on the shock
                RRR=zeros(nb,nb);               
                RRR(:,:)=RR(:,:,l);                     % Take the U(k,k') matrix
                zz=RRR+beta*(v*Py(l,:)')*ones(1,nb);  % This is Bellman equation
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
        bdecis = zeros(nb,ny); 
        cdecis = zeros(nb,ny);
        %idecis = zeros(nb,nstates);
        %indecis = zeros(nb,nstates);
        for i1=1:size(tdecis,1)
            for i2=1:size(tdecis,2)
                bdecis(i1,i2)  = b_grid(tdecis(i1,i2));
                    % Notice that kdecis containts the same information as
                    % tdecis (see above comment), but now it tells you the real
                    % value of i2' (a real number) that is optimal, 
                    % not the position of the grid
                cdecis(i1,i2)  = y_grid(i2)+b_grid(i1)-q(i1)*bdecis(i1,i2);
                %idecis(i1,i2)  = bdecis(i1,i2)- (1-delta)*bgrid(i1);
                %indecis(i1,i2) = idecis(i1,i2) - delta*bgrid(i1);
            end
        end


% ======================================================================
% =========   NEW INTEREST RATE  =======================================
% ======================================================================
bL=bdecis(:,1);
bH=fliplr(bdecis(:,2)')';

mktc=bH+bL;
ssr=sum(mktc.^2);

end
