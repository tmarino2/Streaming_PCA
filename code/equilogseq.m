%% Function generates a sequence close to uniform grid on semilog axis
function [seq,L]=equilogseq(N,numpasses)
L=zeros(1,numpasses+1);
seq=[];
for i=1:numpasses
    tseq=logseq(N*i);
    seq=sort(union(seq,tseq),'ascend');
    L(i+1)=length(seq);
end
end
