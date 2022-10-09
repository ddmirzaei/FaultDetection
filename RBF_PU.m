function s = RBF_PU(X, Y, IndXY, IndYY, fX, rho, op)
% This function implements the RBF-PU method with the constant-generated weight
% Inputs:
%  X: set point of size (N x d)
%  Y: patch centers of size (Nc x d)
%  fX: f values at X of size (N x 1)
%  rho: vector of patch radiuses
%  op: operator
% Output:
%  s: approximate vector (or cell) of size Nx1
%%
Nc = size(Y,1); 
numop = length(op);
if numop ==1
    op = {op};
end
ScaleOrd = ScalingOrder(op);
s = cell(1,numop);
for k = 1:numop
  s{k} = zeros(size(X,1),1);
end
for j = 1:Nc
  yj = Y(j,:);               % j-th patch center   
  Jj = IndXY(j,:);           % indices of points in the j-th patch
  Xj = X(Jj,:);              % points in the j-th patch
  D = DistMat(Xj,Y(IndYY{j},:)); 
  [Dmin,Dind] = min(D,[],2); % find minimum distances
  % indices of points in which yj is their closest center
  inds = (Dind == 1);  
  indx = Jj(inds);           % going from local indices to global ones
  if ~isempty(inds)
    Xhat = Xj(inds,:);       % points in which yj is their closest center
    fXj = fX(Jj);            % function values in the j-th patch
    Lmat = LagrangeMat(Xhat,Xj,yj,rho(j),op,ScaleOrd); % RBF Lagrange matrix
    for k=1:numop
        sj = Lmat{k}*fXj;
        s{k}(indx) = sj;     % putting local approximations sj to global vector s
    end
  end
end
if numop == 1
    s = s{1};
end
