function XNarrow = PCR_PU(X, Y, PatchRad)
% This function combines PCR and PU for narrowing a cloud of points
% Inputs:
%  X: Cloudes of points tat we want to narrow them.
%  Y: Patches center.
%  PatchRad: Radius of patches
% Output:
%  XNarrow: Narrowed points.
%%
XNarrow = zeros(size(X,1),2);
IndXY = PointsInPatch(X,Y,PatchRad);
IndYY = PointsInPatch(Y,Y,1.5*PatchRad);   
for j = 1:size(Y,1)
    ind = IndXY{j}; indc = IndYY{j};
    Xj = X(ind,:);
    D = DistMat(Xj,Y(indc,:));
    [Dmin,Dind] = min(D,[],2);       
    inds = ((Dind == 1));            
    indx = ind(inds);                
    if ~isempty(indx)
        mu = mean(Xj);
        Xj = Xj-mu;                    
        [U,S,V] = svd(Xj','econ');     
        Xpcr = U(:,1)*S(1,1)*V(inds,1)' + mu';
        XNarrow(indx,:) = Xpcr';        
    end
end
