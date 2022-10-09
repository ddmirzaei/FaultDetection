function X = ScatPoints2D(a,b,PointType)
% X: scattered points of size (N x 2).
% a,b: Omega = [a,b]^2.
% PointType: the type of point distribution

switch (PointType)
    case('Halton10000')              
        H = haltonset(2,'Skip',1e4,'Leap',1e2);
        X = net(H,10000); X = (b-a)*X+a;
    case ('RandPoints10000')       
        X0 = load('RandPoints10000'); X = X0.X;        
    case ('RandPoints5000')        
        X0 = load('RandPoints5000'); X = X0.X;
    case 'VarDensityRandPoints'
        X0 = load('VarDensityRandPoints'); X = X0.X;        
    case ('ToyExamplePoints')    
        H = haltonset(2,'Skip',1e4,'Leap',1e2);
        X = net(H,30000); X = (b-a)*X+a;
end
