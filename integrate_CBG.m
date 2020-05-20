function y = integrate_CBG(P,M,U)

x     = M.x;
n     = numel(x);
N     = U.N;
dt    = U.dt;
time  = 0:dt:dt*(N-1);
input = U.input; 

% the output non-linearity
g     = spm_funcheck(M.g);

%% get intrinsic and extrinsic delayed interactions

nmm   = [4 3];      % signify that this is MMC and then BGT

% intrinsic delays ms (scaling)
%--------------------------------------------------------------------------
D(1) = 2;                                % ERP connections      
D(2) = 1;                                % CMC connections
D(3) = 4;                                % BGT connections 
D(4) = 1;                                % MMC connections

% extrinsic delay ms
%--------------------------------------------------------------------------
de   = 8;
di   = 1;

for i = 1:n
    P.D(i,i) = P.D(i,i) + log(D(nmm(i)));
end

De  = exp(P.D);
Di  = diag(diag(De));
De  = De - Di;
De  = De*de/1000;
Di  = Di*di/1000;
D   = Di+De;

Di  = round(diag(D)/dt);
De  = round([D(1,2) D(2,1)]/dt);

% Euler update scheme accounting for delays
% 
neg_t = ceil(max(D(:))/dt);
xs    = repmat(spm_vec(x),1,neg_t+1);

for t = 1:numel(time)
    u = input(t,:);
    temp = [];
    for s = 1:size(xs,1)
        x = spm_unvec(xs(:,neg_t+t),x);
        
        if s<=numel(x{1}), i=1; else i=2; end
          
        if i==1
           xd       = setdiff(1:numel(x{1}),s);
           x{1}(xd) = xs(xd,neg_t+t - Di(1));
           x{2}(:)  = xs(numel(x{1})+1:end,neg_t+t  - De(1));
        else
           ss       = mod(s,numel(x{1}));
           xd       = setdiff(1:numel(x{2}),ss);
           x{1}(:)  = xs(1:numel(x{1}),neg_t+t  - De(2));
           x2       = xs(numel(x{1})+1:end,neg_t+t - Di(2));
           x{2}(xd) = x2(xd); 
            
        end  
        
        % update each state separately and store
        f = spm_fx_gen_ash(x,u,P,M);
        temp(:,s) = spm_vec(x) + f*dt;       
    end
    % update over states
    xs(:,neg_t+t+1)  = diag(temp);
    
    % output - implement g(x)
    y(:,t) = g(xs(:,neg_t+t+1),u,P,M);
    
end



