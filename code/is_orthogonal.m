
%% is_orthogonal_set(A) returns 1 if a set of vectors is orthogonal. 
function FLAG=is_orthogonal(A)
[~,n]=size(A);
TOL=1e-10;
FLAG=1;
if(n==1)
    return;
else
    for i=1:n
        for j=1:i-1
            if(abs(A(:,i)'*A(:,j))>TOL)
                FLAG=0;
                return;
            end
        end
    end
end
end
