function u = FindDirection(X)
% This function finds the principal direction of a set of points by PCR 
% Input:
%  X: points 
% Output:
%  u: principal direction
%%
mu = mean(X);
X = X-mu;                 % zero-mean data
[U,~,~]=svd(X','econ');   % thin svd
u = U(:,1)';
