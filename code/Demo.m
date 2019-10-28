addpath(genpath('../../MBMSG'))

%% Init parameters
k=7;

% RST=0; % do not start from scratch, continue on previous run!
RST=1; % start from scratch, do not continue on previous run!

numiters=1; maxiter=Inf;

eta = .5;
batch_const = 10;

methods={'msg','l2rmsg','sgd','vrmsg','l2vrmsg'};
method_names = {'MSG','RMSG','Oja','VR-MSG','VR-RMSG'};


%% Arbitrary Dataset
% d=1000; gap = 0.1; num_points = 4*1e4;
% syn_gen_exp_gap( k,d,gap, num_points ) 
load('../data/syn_k=7_d=1000_gap=1.00e-01_num_points=40000.mat'); dataname='syn_exp_decay'; tau=0;
% load('sMNIST.mat'); dataname='sMNIST'; tau=0;

%% normalize the dataset
N=size(data.training,2);
N=min(N,maxiter);
data.training=data.training(:,1:N);
mu=mean(data.training,2);
data.training=data.training-repmat(mu,1,N);
data.tuning=data.tuning-repmat(mu,1,size(data.tuning,2));
data.testing=data.testing-repmat(mu,1,size(data.testing,2));

% Renormalize the data
[d, N]=size(data.training);
lambda_k = svds(data.training/sqrt(N),k);
data.training = data.training/lambda_k(k);

S = svds(data.training/sqrt(N),k+1).^2;
S = sort(S,'descend');
gap = S(k) - S(k+1); 
lambda = gap/2;
beta=(S(k)+S(k+1))/4;
%% run the algorithms
stochPCA(dataname,methods,data,k,numiters,eta,batch_const,lambda,RST,beta);

%% plot the results
plotobjV(dataname,methods,k,N,numiters,eta,lambda,beta,method_names);
