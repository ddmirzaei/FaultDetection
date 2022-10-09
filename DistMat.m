function DM = DistMat(Xe,X)
% This function computes the distance matrix of points Xe and X
% Inputs:
%   Xe: evaluation points of size (M x d)
%   X: trial points of size (N x d)
% Outputs:
%   DM: distance matrix of size (M x N)
[n,dim] = size(X); m = size(Xe,1); DM = 0;
for d = 1:dim
  DM = DM + (repmat(X(:,d)',m,1) - repmat(Xe(:,d),1,n)).^2;
end
DM = sqrt(DM);
end % end of function

