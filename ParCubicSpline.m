function s = ParCubicSpline(X,N,SmoothPar)
% This function gives the parametric cubic spline interpolation on ordered 
%     points X in a plane using 'csaps' function of Matlab
% Inputs:
%  X: Ordered points
%  N: Number of evalution points
%  SmoothPar: % Smoothing parameter
% Outputs:
%  s: spline values at evalution point
%%
t (1) = 0;                         
for i = 1:size(X,1)-1
t(i+1) = t (i) + sqrt ( (X(i+1,1)-X(i,1))^2 + (X(i+1,2)-X(i,2))^2 );
end
z = linspace(t(1),t(end),N);
sx = csaps (t,X(:,1),SmoothPar,z);   
sy = csaps (t,X(:,2),SmoothPar,z);
s = [sx' sy'];
