
%% incremental updates the EIG of the covariance matrix Cxx given
% a new sample x using an incremental SVD approach
% Usage:
% Inputs:  U,S   - current EIG
%          x     - new sample pair
%          RANK  - maximum EIG-rank desired
% Outputs: U,S   - updated EIG
%
%%
function [U,S]=incremental(U,S,x,RANK)
xproj=U'*x;          % Project data onto left singular subspace (rx1)
xres=x-U*xproj;      % Residual unexplained by left singular subspace (mx1)
xresnorm=norm(xres); % Norm of the residual
if(xresnorm>1e-6)
    xres=xres/xresnorm; % Normalize residual to get unit vector normal to U
end

Q=[S+(xproj*xproj') xresnorm*xproj;...
    xresnorm*xproj' xresnorm^2];% Form the (r+1)x(r+1) matrix for inner SVD
[U2,S2,~]=svd(Q,'econ');%Inner SVD: U2, S2, V2 are all (r+1)x(r+1) matrices
U=[U xres]*U2;       % U is mx(r+1).
S=S2;                % S is (r+1)x(r+1)
if(nargin>3)         % TRUNCATE the SVD to the given RANK
    if(length(diag(S>0))>RANK)%Check if nonzero singular values exceed RANK
        dS=diag(S);  % Singular values
        [~,idx]=sort(dS,'descend');%Sort singularvalues in descending order
        idx=idx(1:RANK);  % Only top RANK singular values, vectors needed
        U=U(:,idx);       % Top left-singular vectors
        S=diag(dS(idx));  % Top singular values
    end % End SVD truncation
end % End-IF
end % End-function
