
%% Performs the "positive" Frobenius-MD update. The iterate M is:
%     M = U * diag( S ) * U'
%
% k - the dimension of the subspace which we seek
% U, S - "nontrivial" eigenvectors and eigenvalues of the iterate
% eta - the step size
% x - the sample vector
% eps - threshold for rank1update and msgproject
%%
function [U,S]=msg(k,U,S,eta,x,eps)

[U,S]=rank1update(U,S,eta,x,eps);

nz=S>0; S=S(nz); U=U(:,nz);
[S,idx]=sort(S,'ascend');
U=U(:,idx);
% S is sorted in ascending order
S=the_projection(S,k);
end
