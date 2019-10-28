
%% GETREF(A) returns pivot columns from reduced echelon form
function U=getref(A)
[m,n]=size(A);
flag = 0;
if(n==2)
    multiple=A(1,2)/A(1,1);
    count=0;
    for i=1:m
        if(A(i,2)/A(i,1) == multiple)
            count=count+1;
        end
    end
    if(count==m)
        U=A(:,1);
        flag=1;
    end
end
if(flag==0)
    [~,pivot_columns]=ref(A);
    for i=1:size(pivot_columns,2)
        U(:,i)=A(:,pivot_columns(1,i));
    end
end
end
