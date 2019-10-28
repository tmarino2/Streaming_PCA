%% Wrapper for running the experiments for VRMSG
% k is number of components to be recovered
% numiters is number of times to repeat the experiment
% dataname is the name of the dataset
% d is dimensionality of data (provide only if generating synthetic data)
% gap is as described in Section 7 of the paper (provide only if generating synthetic data)
% num_points is the number of points in the dataset (provide only if generating synthetic data)

function vrmsgWrapper (k,numiters,dataname,d,gap,num_points)
    addpath(genpath('../../MBMSG'));
    if(strcmp(dataname,'synthetic'))
        syn_gen_exp_gap(k,d,gap, num_points)
        fname = sprintf('../data/syn_k=%d_d=%d_gap=%0.2e_num_points=%d.mat',k,d,gap,num_points);
        load(fname);
    else
        fname = sprintf('../data/%s.mat',dataname);
        load(fname);
    end    
    methods={'msg','l2rmsg','sgd','vrmsg','l2vrmsg'};
    method_names = {'MSG','RMSG','Oja','MB-MSG','MB-RMSG'};

    eta = .5;
    batch_const = 10;
    RST = 1;
    
    N=size(data.training,2);
    data.training=data.training(:,1:N);
    mu=mean(data.training,2);
    data.training=data.training-repmat(mu,1,N);
    data.tuning=data.tuning-repmat(mu,1,size(data.tuning,2));
    data.testing=data.testing-repmat(mu,1,size(data.testing,2));
    lambda_k = svds(data.training/sqrt(N),k);
    data.training = data.training/lambda_k(k);
    
    
    S = svds(data.training/sqrt(N),k+1).^2;
    S = sort(S,'descend');
    gap = S(k) - S(k+1); 
    beta = (S(k) + S(k+1))/4;
    lambda = gap/2;
    
    stochPCA(dataname,methods,data,k,numiters,eta,batch_const,lambda,RST,beta);
    
    plotobjV(dataname,methods,k,N,numiters,eta,lambda,beta,method_names);
end

