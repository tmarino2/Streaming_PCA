
%% Gram_Schmidt(A) produces an orthonormal basis for the subspace spanned by 
% the vectors A=[a1,a2,...,an] using Gram-Schmidt's orthogonalization
% procedure
function [U]=Gram_Schmidt(A)

[m,n]=size(A);
if(norm(A-zeros(m,n))<1e-10)
    error('Zero Vector Basis: no basis exist');
elseif(n==1)
    U=A(1:m,1)/norm(A(1:m,1));
else
    if(is_orthonormal(A)==1)
        U=A;
        return;
    end
    if(rank(A)~=n)
        A=getref(A);
    end
    [m,n]=size(A);
    U=A(1:m,1)/norm(A(1:m,1));
    for i=2:n
        u=A(1:m,i);
        v=u;
        for j=1:(i - 1)
            v=v-(u'*U(1:m,j))*U(1:m,j);
        end
        U(:,i)=v/norm(v);
    end
end
end
