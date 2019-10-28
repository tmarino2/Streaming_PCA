%% Performs the inefficient MB-MSG (MB-RMSG) update:
% Use only if you want to verify the projection returns a rank k matrix!
% k - the dimension of the subspace which we seek
% U, S - "nontrivial" eigenvectors and eigenvalues of the iterate
% eta - the step size
% X - matrix with columns the samples for the t-th iteration
% eps - threshold for msgproject (discard directions with eigenvalue<eps)
%%
function [U,S]=vrmsg(k,U,S,eta,X,eps)
    n = size(X,2);
    X_temp = sqrt(eta/n)*X;
    [P,~,~] = svds(X_temp - U*(U'*X_temp));
    R_X = P'*(X_temp - U*(U'*X_temp));
    UU = [U P];
    kk = size(P,2);
    Upd = [U'*X_temp; R_X];
    K = [diag(S), zeros(k,kk); zeros(kk,k) zeros(kk,kk)] + Upd*Upd';
    [Rot, S] = eig(K);
    U = UU*Rot;
    S = diag(S);
%     C_t = 1/n*(X*X');
%     P_t12 = U*diag(S)*U' + eta*C_t;
%     [U,S] = eig(P_t12);
%     S = diag(S);

    % S is sorted in ascending order  
    [S,idx]=sort(S,'ascend');
    U=U(:,idx);    
    [U,S]=projection(U,S,k);
    U = U(:,S>eps);
    S = S(S>eps);

%     nz=S>0; S=S(nz); U=U(:,nz);
%     [S,idx]=sort(S,'ascend');
%     U=U(:,idx);
%     % S is sorted in ascending order
%     S=the_projection(S,k);
%     
%     U = U(:,S>eps);
%     S = S(S>eps);    
end
