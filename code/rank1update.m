
%% Performs the rank-1 update:
%     M + eta * sample * sample'
% where M is represented as a partial eigendecomposition:
%     M = U * diag(S) * U'
%
% U, S - "nontrivial" eigenvectors and eigenvalues of M
% eta - the scale of the rank-1 update
% sample - the vector of the rank-1 update
% eps - threshold for adding a new "nontrivial" dimension
%%
function [U,S]=rank1update(U,S,eta,x,eps )
[d,k]=size(U);
xhat=U'*x;
res=x-U*xhat;
resnorm=norm(res);
if ((resnorm < eps) || (k >= d))
    dS=diag(S)+eta*(xhat*xhat');
    dS=0.5*(dS+dS'); %**make sure matlab knows the matrix is real symmetric
    
else
    U=[U,res/resnorm ];
    dS=[diag(S)+eta*(xhat*xhat'), (eta*resnorm)*xhat;...
        (eta*resnorm)*xhat', (eta*resnorm*resnorm)];
    dS=0.5*(dS+dS'); %**make sure matlab knows the matrix is real symmetric
end
[Utilde,newS]=eig(dS);
U=U*Utilde;
S=diag(newS)';
U=real(U);
S=real(S);
end
