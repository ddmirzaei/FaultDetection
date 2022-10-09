function K = KerMat(X,Y,op)
% This function computes the kernel matrix for different operator 
%  on two point sets X and Y
% Inputs:
%   X: test points of size (m x dim)
%   Y: trial points of size (n x dim)
%   op: operator: '0' for identity, 'x' for first derivative respect to x,
%       etc ...
% OutputsL    
%   K: kernel matrix of size (m x n)

dim = size(X,2); % dimension
s = DistMatSqH(X,Y);
switch (op)
    case('1')
        K = Frbf(s,0);
    case('x')
        x = DiffMat(X(:,1),Y(:,1));
        K = x.*Frbf(s,1);
    case('y')
        y = DiffMat(X(:,2),Y(:,2));
        K = y.*Frbf(s,1);
    case('z')
        z = DiffMat(X(:,3),Y(:,3));
        K = z.*Frbf(s,1);
    case('xx')
        x = DiffMat(X(:,1),Y(:,1));
        K = Frbf(s,1)+x.^2.*Frbf(s,2);
    case('yy')
        y = DiffMat(X(:,2),Y(:,2));
        K = Frbf(s,1)+y.^2.*Frbf(s,2);
    case('zz')
        z = DiffMat(X(:,3),Y(:,3));
        K = Frbf(s,1)+z.^2.*Frbf(s,2);
    case ('L')   %\Delta
        K = dim*Frbf(s,1)+2*s.*Frbf(s,2);
    otherwise
        error('this type of KerMat operator (char argument) is not implemented')
end