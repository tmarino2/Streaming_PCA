%% Performs the fast MB-MSG update.:
% This function is directly implemented in the main code stochPCA.m
% % k - the dimension of the subspace which we seek
% U, S - "nontrivial" eigenvectors and eigenvalues of the iterate
% eta - the step size
% X - matrix with columns the samples for the t-th iteration
%%
function [U,S]=fast_vrmsg(k,U,S,eta,X)
    n = size(X,2);
    X_temp = sqrt(eta/n)*X;
    [U,~,~] = svds([U,X_temp],k);  
end
