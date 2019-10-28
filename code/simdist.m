
%% SIMDIST draws samples from the orthogonal distribution described
% in homework problem 5
%
% Usage:    [X] = simdist(d,n,tau)
% Input:    d    - dimensionality (32 in homework)
%           n    - number of samples to be drawn
%           tau  - spectrum roll-off parameter
% Output:   X    - dxn matrix containing samples drawn as columns
%
% Example usage: X=simdist(32,100000,1.1);
%
% Version History
% Ver 0.1 Created Oct 5th, 2014
% Author: Raman Arora
%%

function [X,sigma]=simdist(d,n,tau)
sigma=(tau).^(-(1:d));
sigma=sigma./sum(sigma);
I=eye(d);
S=sign((rand(1,n)>0.5)-0.5);
X=repmat(S,d,1).*I(:,sampA(sigma,n));
DEBUG=0;
if(DEBUG)
    mx=mean(X,2); %#ok<UNRCH>
    fprintf('Norm of the mean vector: %g\n',norm(mx));
    cx=cov(X');
    fprintf('Deviation from the desired spectrum: %g\n',norm(diag(cx)-sigma'));
    figure; clf; plot(sigma,'Linewidth',2); hold on; plot(diag(cx),'r','Linewidth',2);
    legend('Desired spectrum','Realized spectrum'); set(gca,'FontSize',16);
    grid; axis([0 d -0.2 0.25]);
end
end
