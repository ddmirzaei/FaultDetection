function Ind = PointsInPatch(X,Y,rad)
% X: a set of points 
% Y: patch's centers
% rad: radius of patches.
% Ind: indices of points in each patch
T = KDTreeSearcher(X);   
Ind = rangesearch(T,Y,rad);
end % end of function
