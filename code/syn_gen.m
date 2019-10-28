% clear all
% X=randn(100,40000);
% [U,~,V]=svd(X,'econ');
% S=diag(1.2.^(0:-.5:-49.5));
% X=U*S*V';
% m=mean(X,2);
% X=X-repmat(m,[1,40000]);
% s=0;
% for i=1:40000
%     s=max(s,norm(X(:,i)));
% end
% X=X/s;
% data.training=X(:,1:20000);
% data.tuning=X(:,20001:30000);
% data.testing=X(:,30001:40000);
% save synthetic.mat data
% 


%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%% big gap at 2,4,8

% clear all
% X=rand(100,10000);
% [u s v]=svd(X,'econ');
% s=sqrt([1,.98,.68,.67,.37,.36,.35,.34,.04]);
% s=[s,s(9):-.00219:0];
% s(9)=[];
% X=u*diag(s)*v'*100;
% data.training=X(:,1:6000);
% data.tuning=X(:,6001:8000);
% data.testing=X(:,8001:10000);
% save syn.mat data

%% exp decay 

nn=4*1e4;
% X=rand(100,nn);
% [uu ss vv]=svd(X,'econ');
% ss=2.^(0:-.1:-9.9);
% X=uu*diag(ss)*vv'*sqrt(nn);

d=1000;
mu = zeros(1,d);
step = .4;
Sigma = 2.^(0:-step:-(d-1)*step);
%Sigma = 0.1*ones(1,d);
Sigma(1:7) = 1;
X = mvnrnd(mu,Sigma,nn);
X = X';
data.training=X(:,1:.8*nn);
data.tuning=X(:,.8*nn+1:.9*nn);
data.testing=X(:,.9*nn+1:nn);
save ../data/syn_exp_decay.mat data
