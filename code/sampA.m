
%% Function sampA draws multiple samples from a specified pmf A
%
% Usage:    [kN] = sampA(A,N)
% Input:    A    - Specified pmf
%           N    - Number of sampled to be drawn
% Output:   kN   - Samples drawn from the pmf
%
% Version History
% Ver 0.1 Created Oct 5th, 2014
% Author: Raman Arora

function [kN] = sampA(A,N)
kN = zeros(1,N);
lenA = length(A)+1;
cumA = zeros(1,lenA);
cumA(2:end) = cumsum(A);
krand = rand(1,N);
for i = 2:length(cumA)
    relem = intersect(find(krand > cumA(i-1)),find(krand <= cumA(i)));
    kN(relem) = i-1;
end
end
