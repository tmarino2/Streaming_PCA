%% Checks if the sufficient condition for the projection is satisfied on the current sample
% k is the number of components
% U is the matrix of eigenvectors of P_t
% X is the current sample
% gap_const is 2 or 4 depending on vrmsg or vrrmsg
function [sat]=suff_cond(k,U,X,gap_const)  
    Sigma = svds(X,k+1).^2;
    Sigma = sort(Sigma,'descend');
    gap = Sigma(k) - Sigma(k+1);
%     sat = (sum(Sigma(1:k)) - trace((U'*X)*(X'*U)) > gap/gap_const);
    lambda_k = svds(U'*X,k).^2;
    lambda_k = sort(lambda_k,'descend');
    lambda_k = lambda_k(k);
    if(gap_const == 2)
        sat = ( trace((U'*X)*(X'*U)) + lambda_k - sum(Sigma) < 0);
    else
        sat = ( trace((U'*X)*(X'*U)) + lambda_k - sum(Sigma) - gap/2< 0);
    end
end
