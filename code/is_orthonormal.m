
%% is_orthonormal(A) returns 1 if a set of vectors is orthonormal. 
function FLAG = is_orthonormal(A)
[m,n]=size(A);
TOL = 1e-10;
FLAG=1; 
if(norm(A-zeros(m,n))<1e-10)
    error('Zero Vector Basis: no basis exist');
elseif(n==1)
    if(abs(norm(A)-1)>TOL)
        FLAG=0;
        return;
    end
else
    for i=1:n
        if(abs(norm(A(:,i))-1)>TOL)
            FLAG=0;
            return;
        end
    end
	FLAG=is_orthogonal(A);
end
end
