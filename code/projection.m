%% Efficient implementation of the projection \Pi in the paper.
% U, S - "nontrivial" eigenvectors and eigenvalues of the iterate
% k - the dimension of the subspace which we seek
function [ U, S ] = projection( U, S, k )
%Assume S is a vector sorted in ascending order
d = length(S);

% check if capping without shift is OK
SS=S;
SS(S>1)=1; SS(S<0)=0;
if sum(SS)<=k
    S=SS;
    return;
end

% all negative entries are going to be capped to 0 anyway
S(S<0) = 0;

SS = S;
i=1; j=1;
while 1-SS(j) > 0 && j < d
    j = j+1;
end

while j<d-1
    s_ij = 1-SS(j+1);
    while SS(i) + s_ij <= 0 && i <= j+1
        i = i+1;
    end
    if i>j+1
        j = j+1;
    end
    p_ij = sum(max(0,min(1,SS+s_ij)));
    if p_ij > k
        j = j+1;
        continue
    end
    if p_ij == k
        %S = SS;
        S = max(0,min(1,SS+s_ij));
        return;
    end
    %find i
    i = 1;
    while i<=j
        s_ij = (k-(d-j)-sum(S(i:j)))/(j-i+1);
        if s_ij <0 && SS(i) + s_ij >= 0 && ((i > 1 && S(i-1) + s_ij <= 0) || i==1)...
                && SS(j) + s_ij <= 1 && S(j+1) + s_ij >= 1
            S = max(0,min(1,SS+s_ij));
            return;
        end
        i = i+1;
    end
end
for j = [d-1,d]
i = 1;
    while i<=j
        s_ij = (k-(d-j)-sum(S(i:j)))/(j-i+1);
        if s_ij < 0 && SS(i) + s_ij >= 0 && ((i > 1 && S(i-1) + s_ij <= 0) || i==1)...
            && SS(j) + s_ij <= 1 && ((j<d && S(j+1) + s_ij >= 1) || j==d)
            S = max(0,min(1,SS+s_ij));
            return;
        end
        i = i+1;
    end
end
S(1:d-k) = 0;
S(d-k+1:d) = 1;
end

