%% SGD is the stochastic gradient descent algorithm for PCA.
% It updates the parameter basis matrix U based on the estimate
% of the gradient xx' evaluated on the new sample x
% 
% Usage:
% Inputs:  U     - current basis 
%          x     - new sample
%          etax  - step-size
% Outputs: U	 - updated basis (not orthogonalized)
%
%%
function [U]=sgd(U,x,etax)

U=U+etax*(x*(x'*U)); % Projection step needed infrequently

end % End-function

