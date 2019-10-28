function [ U,S ] = msgcapping( U,S,K )
%MSGCAPPING Summary of this function goes here
%   Detailed explanation goes here
if length(S)>K
    [S,srt]=sort(S,'descend');
    S=S(1:K);
    U=U(:,srt(1:K));
end
end

