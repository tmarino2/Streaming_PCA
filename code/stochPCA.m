
%%STOCHPCA(dataname,methods,data,k,K,iter_end,eta,lambda,RST,beta,iter_start)
% Runs various stochastic PCA algs
% The inputs are as follows:
%  dataname  - 'orthogonal' or 'mnist'
%  methods   - 'batch', 'incremental', 'msg', 'l1rmsg', 'l2rmsg', 'l21rmsg'
%  data      - data.training, data.tuning, data.testing
%   k        - a positive integer, denotes the desired RANK
%   K        - capping rank
%   iter_end - end of iterations, iter_start ... iter_end
%   eta      - step size multiplier. The sequence of step size used by
%             various algorithms is eta/t or eta/sqrt(t), where t is the 
%             iteration index.
%   lambda   - l2 parameter
%   RST      - 0 -> continue, 1 -> restart
%   beta     - l1 parameter
%   itr_st   - start of iterations, iter_start ... iter_end
%
%  The output containing the objective on the dev set,population objective,
%  and the singular value decomposition (U,S,V) is written to a
%  MATFILE in ../page/profile/PCA/METHOD/DATANAME, where METHOD is one of
%  'sgd', 'incremental', 'msg', 'batch', or 'truth'. For 'orthogonal'
%  distribution, the filename follows the format:
%  'method[rank=1,tau=1,iter=1].mat'
%%
function stochPCA(dataname,methods,data,k,iter_end,eta,batch_const,lambda,RST,beta,itr_st)


%% Default is starting at iter_start=1;
if(nargin<11)
    itr_st=1;
end

%% Default is beta=.5
if(nargin<10)
    beta=0.5;
end


%% Default is not to restart
if(nargin<9)
    RST=0;
end


%% Default is lambda = 1e-3
if(nargin<8)
    lambda=1e-3;
end


%% Default is no step size = 1
if(nargin<7)
    eta=1;
end

%% Default is just one run
if(nargin<6)
    iter_end=1;
end

%% Computing Cov matrices
CXX=(1/(size(data.testing,2)-1))*...
    (data.testing)*...
    (data.testing)';

CXXtune=(1/(size(data.tuning,2)-1))*...
    (data.tuning)*...
    (data.tuning)';


%% Simulation parameters
[d, N]=size(data.training);
% % Renormalize the data
% lambda_k = svds(data.training/sqrt(N),k);
% data.training = data.training/lambda_k(k);

%% Check if all the runs are done
flag=1;
for method=methods
    %% IF reading from or writing to a report file
    if(sum(strcmp(method,{'batch','incremental'})))
        reportfile=sprintf(['../report/reportPCA[data=%s,method=%s,',...
            'rank=%d,numiter=%d].mat'],...
            dataname,method{1},k,iter_end);
    elseif(sum(strcmp(method,{'l2rmsg'})))
        reportfile=sprintf(['../report/reportPCA[data=%s,method=%s,',...
            'rank=%d,eta=%f,lambda=%f,numiter=%d].mat'],...
            dataname,method{1},k,eta,lambda,iter_end);
    elseif(sum(strcmp(method,{'l1rmsg'})))
        reportfile=sprintf(['../report/reportPCA[data=%s,method=%s,',...
            'rank=%d,eta=%f,beta=%f,numiter=%d].mat'],...
            dataname,method{1},k,eta,beta,iter_end);
    elseif(sum(strcmp(method,{'l21rmsg'})))
        reportfile=sprintf(['../report/reportPCA[data=%s,method=%s,',...
            'rank=%d,eta=%f,lambda=%f,beta=%f,numiter=%d].mat'],...
            dataname,method{1},k,eta,lambda,beta,iter_end);
    else
        reportfile=sprintf(['../report/reportPCA[data=%s,method=%s,',...
            'rank=%d,eta=%f,numiter=%d].mat'],...
            dataname,method{1},k,eta,iter_end);
    end
    flag=flag && (exist(reportfile,'file')) && ~RST;
end

if(~flag)
    for ITER=itr_st:iter_end
        
        %% Shuffle data
        rng(ITER);
        rp=randperm(N);
        data.training=data.training(:,rp);
        
        %% Sequence close to a uniform grid on semilog axis
        [seq,L]=equilogseq(N,1);
        
        for method=methods
            
            %% Display the run
            fprintf('Starting run: (%s,%s,%s,%d,%d)\n',...
                dataname,'PCA',method{1},k,ITER);
            
            %% Set PAGE directories
            pagepath=sprintf('../page/profile/PCA/%s/%s/',...
                method{1},dataname);
            
            if(sum(strcmp(method,{'batch','incremental'})))
                pageprefix=@(method,rank,iter)[pagepath,...
                    sprintf('%s[rank=%d,iter=%d].mat',...
                    method,rank,iter)];
            elseif(sum(strcmp(method,{'l2rmsg'})))
                pageprefix=@(method,rank,iter)[pagepath,...
                    sprintf('%s[rank=%d,eta=%f,lambda=%f,iter=%d].mat',...
                    method,rank,eta,lambda,iter)];
            elseif(sum(strcmp(method,{'l1rmsg'})))
                pageprefix=@(method,rank,iter)[pagepath,...
                    sprintf('%s[rank=%d,eta=%f,beta=%f,iter=%d].mat',...
                    method,rank,eta,beta,iter)];
            elseif(sum(strcmp(method,{'l21rmsg'})))
                pageprefix=@(method,rank,iter)[pagepath,...
                    sprintf('%s[rank=%d,eta=%f,lambda=%f,beta=%f,iter=%d].mat',...
                    method,rank,eta,lambda,beta,iter)];
            else
                pageprefix=@(method,rank,iter)[pagepath,...
                    sprintf('%s[rank=%d,eta=%f,iter=%d].mat',...
                    method,rank,eta,iter)];
            end
            % Check if the PAGE directory is structured properly
            if(~exist(pagepath,'dir'))
                % If not create the desired directory structure
                flag=createpath(pagepath);
                % If the directory structure could not be created
                if(~flag)
                    % Display error message and quit
                    error('Could not create path for result files');
                end
            end
            
            %% Output filename
            fname=pageprefix(method{1},k,ITER);
            if(~exist(fname,'file') || RST)
                
                CXX=(1/(size(data.testing,2)-1))*...
                    (data.testing)*...
                    (data.testing)';
                
                CXXtune=(1/(size(data.tuning,2)-1))*...
                    (data.tuning)*...
                    (data.tuning)';
                
                %% Initialize the basis for sgd and incremental PLS
                gap = 1;
                S=0;
                if(strcmp(method{1},'sgd'))
                    U=orth(randn(d,k));
                elseif(sum(strcmp(method{1},{'incremental','msg','l2rmsg','l1rmsg','l21rmsg'}))) %for msg and l2rmsg start at zero, project wrt trace<=k
                    U=zeros(d,k);
                    gap = 2*lambda;
                    S=zeros(k,1);
                elseif(strcmp(method{1},{'vrmsg'}))
                    gap = 2*lambda;
                    mini_batch = ceil(batch_const*(k^2)/(gap^(3)));
                    U=orth(randn(d,k));
                    S=ones(k,1);
                    indx_set = randi(N,1,min(mini_batch,N));
                    X_init_test = data.training(:,indx_set);
                    indx_set = randi(N,1,min(mini_batch,N));
                    samples_t = k;
                    while(samples_t < mini_batch && suff_cond(k,U,X_init_test/sqrt(mini_batch),2))
                        samples_t = min(samples_t*2,mini_batch);
                        [U,~,~] = svds(data.training(:,indx_set(1:samples_t)),k);
                    end
                    
                elseif(strcmp(method{1},{'l2vrmsg'}))          
                    gap = 2*lambda;
                    mini_batch = ceil(batch_const*(sqrt(k))/(gap^(5)));
                    U=orth(randn(d,k));
                    S=ones(k,1);
                    indx_set = randi(N,1,min(mini_batch,N));
                    X_init_test = data.training(:,indx_set);
                    indx_set = randi(N,1,min(mini_batch,N));
                    samples_t = k;
                    while(samples_t < mini_batch && suff_cond(k,U,X_init_test/sqrt(mini_batch),2))
                        samples_t = min(samples_t*2,mini_batch);
                        [U,~,~] = svds(data.training(:,indx_set(1:samples_t)),k);
                    end               
                end
                
                %% Initialize objective value
                rk=zeros(L(2),1);
                objV=zeros(L(2),1);
                objVtune=zeros(L(2),1);
                runtime=zeros(L(2),1);
                mb_size = zeros(L(2),1);
                
                %% Check if we can start from a previous run
                initsamp=1;
                
                %% Loop over data
                for iter=L(1)+1:L(2)
                    fprintf('Sequence number %d...\t',seq(iter));
                    switch(method{1})
                        case 'batch'
                            %% BATCH PCA
                            isamp=seq(iter);
                            if(isamp==1)
                                continue;
                            end
                            tcounter=tic;
                            Ctrain=(1/(isamp-1))*...
                                ((data.training(:,1:isamp))*...
                                (data.training(:,1:isamp))');
                            [U,S,~]=svds(Ctrain,k);
                            U=U(:,1:k);
                            runtime(iter)=toc(tcounter);
                            rk(iter)=size(S,1);
                            
                        case 'sgd'
                            %% Stochastic gradient descent
                            for isamp=initsamp:seq(iter)
                                modisamp=1+mod(isamp-1,N);
                                etax=eta/(modisamp);
                                tcounter=tic;
                                U=sgd(U, data.training(:,modisamp),etax);
                                U=Gram_Schmidt(U);
                                runtime(iter)=toc(tcounter);
                            end
                            rk(iter)=size(U,2);
                            
                        case 'incremental'
                            %% BRAND's method
                            for isamp=initsamp:seq(iter)
                                modisamp=1+mod(isamp-1,N);
                                tcounter=tic;
                                [U,S]=rank1update(U,S,1,data.training(:,modisamp),1e-6);
                                [U,S]=msgcapping(U,S,k);
                                runtime(iter)=toc(tcounter);
                            end
                            rk(iter)=size(U,2);
                            
                        case {'msg'}
                            %% Matrix Stochastic Gradient
                            for isamp=initsamp:seq(iter)
                                modisamp=1+mod(isamp-1,N);
                                etax=eta/sqrt(modisamp);
                                tcounter=tic;
                                [U,S]=msg(k,U,S,etax,...
                                    data.training(:,modisamp),1e-6);
                                runtime(iter)=toc(tcounter);
                            end
                            fprintf('rank of iterate: %d ',length(find(S)));
                            rk(iter)=length(find(S));
                            
                        case {'vrmsg'}
                            %% Variance Reduced MSG
                            for isamp = initsamp:seq(iter)
                                modisamp=1+mod(isamp-1,N);
                                mini_batch = max(ceil(batch_const*(k+1)/gap^2),k+2);
                                indx_set = randi(N,1,mini_batch);
                                etax = eta/sqrt(modisamp + (k^6)/(gap^2));
                                samples_t = k+1;
                                while(samples_t < mini_batch &&...
                                        suff_cond(k,U,data.training(:,indx_set(1:samples_t))/sqrt(samples_t),2))
                                    samples_t = min(samples_t*2,mini_batch);
                                end  
                                X_temp = data.training(:,indx_set(1:samples_t));
                                X_temp = sqrt(etax/samples_t)*X_temp;
                                opts.tol=1e-5;
                                opts.maxit = 150;
                                opts.initY = U;
                                tcounter=tic;
                                [~,~,U,~] = lmsvd([U,X_temp]',k,opts);
                                runtime(iter) = toc(tcounter);
                            end
                            fprintf('rank of iterate: %d ',length(find(S)));
                            fprintf('minibatch size: %d ',samples_t);
                            rk(iter)=length(find(S));
                            mb_size(iter) = samples_t;
                            
                        case {'l2rmsg'}
                            %% l2 Matrix Stochastic Gradient
                            for isamp=initsamp:seq(iter)
                                modisamp=1+mod(isamp-1,N);
                                etax=eta/(lambda*modisamp);
                                tcounter=tic;
                                [U,S]=msg(k,U,(1-etax*lambda)*S,etax,...
                                    data.training(:,modisamp),1e-6);
                                runtime(iter)=toc(tcounter);
                            end
                            fprintf('rank of iterate: %d ',length(find(S)));
                            rk(iter)=length(find(S));
                            
                        case{'l2vrmsg'}
                            for isamp=initsamp:seq(iter)
                                modisamp=1+mod(isamp-1,N);
                                mini_batch = max(ceil(batch_const*(k+1)/gap^2),k+2);
                                etax=eta/(lambda*(modisamp+1/gap^3));
                                indx_set = randi(N,1,mini_batch);                          
                                samples_t = k+1;
                                while(samples_t < mini_batch &&...
                                        suff_cond(k,U,data.training(:,indx_set(1:samples_t))/sqrt(samples_t),2))
                                    samples_t = min(samples_t*2,mini_batch);
                                end                             
                                X_temp = data.training(:,indx_set(1:samples_t));
                                X_temp = sqrt(etax/samples_t)*X_temp;
                                opts.tol=1e-5;
                                opts.maxit = 150;
                                opts.initY = U;
                                tcounter=tic;
                                [~,~,U,~] = lmsvd([U,X_temp]',k,opts);
                                runtime(iter)=toc(tcounter);
                            end
                            fprintf('rank of iterate: %d ',length(find(S)));
                            fprintf('minibatch size: %d ',samples_t);
                            rk(iter)=length(find(S));
                            mb_size(iter) = samples_t;
                            
                        case {'l1rmsg'}
                            %% l1 Matrix Stochastic Gradient
                            for isamp=initsamp:seq(iter)
                                modisamp=1+mod(isamp-1,N);
                                etax=eta/sqrt(modisamp);
                                tcounter=tic;
                                [U,S]=l1rmsg(U,S,k,etax,...
                                    data.training(:,modisamp),1e-6,beta);
                                runtime(iter)=toc(tcounter);
                            end
                            rk(iter)=length(find(S));
                            
                        case {'l21rmsg'}
                            %% l2+l1 Matrix Stochastic Gradient (elastic net)
                            for isamp=initsamp:seq(iter)
                                modisamp=1+mod(isamp-1,N);
                                etax=eta/(lambda*modisamp);
                                tcounter=tic;
                                [U,S]=l1rmsg(U,(1-etax*lambda)*S,k,etax,...
                                    data.training(:,modisamp),1e-6,beta);
                                runtime(iter)=toc(tcounter);
                            end
                            rk(iter)=length(find(S));
                                                        
                    end
                    
                    initsamp=seq(iter)+1;
                    keff=min(k,size(U,2));
                    if(strcmp(method{1},'sgd'))
                        UU=Gram_Schmidt(U);
                    elseif(sum(strcmp(method{1},{'msg','l2rmsg','l1rmsg', 'l21rmsg','vrmsg','l2vrmsg'})))
                        UU=pca_solution_original(keff,U,S);
                    else
                        UU=U;
                    end
                    
                    if any(UU~=0)
                        UU=Gram_Schmidt(UU);
                    end
                    
                    objVtune(iter)=trace(UU'*CXXtune*UU);
                    objV(iter)=trace(UU'*CXX*UU);
                    fprintf('\t%d\t %g\n',keff,objVtune(iter));
                end
                save(fname,'runtime','U','S','objV','objVtune','seq','rk','mb_size');
            end
        end
        
    end
end


%% Set PAGE directories
method={'truth'};
pagepath=sprintf('../page/profile/PCA/%s/%s/',method{1},dataname);
% Check if the PAGE directory is structured properly
if(~exist(pagepath,'dir'))
    % If not create the desired directory structure
    flag=createpath(pagepath);
    % If the directory structure could not be created
    if(~flag)
        % Display error message and quit
        error('Could not create path for result files');
    end
end
fname=[pagepath,sprintf('truth[rank=%d].mat',k)];
if(~exist(fname,'file'))
    
    [EigVecs,EigVals]=eig(CXX);
    [EigVals, idx]=sort(diag(EigVals),'descend');
    objV=sum(EigVals(1:k)); %#ok<NASGU>
    EigVecs=EigVecs(:,idx(1:k));
    objVtune=trace(EigVecs'*CXXtune*EigVecs); %#ok<NASGU>
    rk=k; %#ok<NASGU>
    save(fname,'objV','objVtune','rk');
    
    
end

end
