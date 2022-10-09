function Y = PatchCenters(PointType,a,b,hcov)
% This function produces patch centers on square [a,b]^2 for all examples
% Inputs:
%  PointType: Type of points 
%  a, b: Determine the domain [a,b]^2 
%  hcov: spacing distance between centers in Y
% Output:
% Y: patch centers 

switch PointType
    case 'VarDensityRandPoints'
        Y0 = load('VarDensityPatches'); Y = Y0.Y;
    otherwise
        y = a:hcov:b;
        [xx,yy] = meshgrid(y,y);
        y1=xx(:); y2=yy(:);
        Y = [y1 y2];
end
