%% Example 1
clear 
%% Function f
f = piecewise(@(x,y) x<=0.7+0.2*sin(2*pi*y),@(x,y) abs(x-0.4-0.1*sin(2*pi*y)),...   
              @(x,y) x>0.7+0.2*sin(2*pi*y),@(x,y) abs(x-0.4-0.1*sin(2*pi*y))-0.2);   
%% Initial point set X and function values
PointType = 'RandPoints10000'; %'RandPoints10000'; 'VarDensityRandPoints'
X = ScatPoints2D(0,1,PointType);  
fX = f(X(:,1),X(:,2));                    
%% Parameter selection
N = size(X,1); % number of initial points 
h = 1/sqrt(N); % approximate fill-distance 
CH = 6;        % constant parameter
H1 = CH*h; H2 = H1; H3 = 0.5*H1; H4 = 2*H1; H5 = H1; % H parameters  
Hpar.H1 = H1; Hpar.H2 = H2; Hpar.H3 = H3; Hpar.H4 = H4; Hpar.H5 = H5; 
CL = 1/2;          % Gradient fault parameter
CG = 1;            % Ordinary fault parameter
CM = 1/4;          % Mean value parameter 
eta = 4*h;         % Parameter for removing false gradient fault points 
Cpar.L = CL; Cpar.G = CG; Cpar.M = CM; Cpar.eta = eta;
SmoothPar = 0.9999;% Smoothing parameter for the cubic smoothing spline.
C_cov = 2.5;
PUpar.hcov = C_cov*h;    % spacing distance between patch centers 
PUpar.ovlp = 1.5;        % overlapping constant (uses in the narrowing step)
%%
Y = PatchCenters(PointType,0,1,PUpar.hcov); % patch centers 
%%
disp('Detection algorithm starts ...')
[NarrowPoints,XnearFault,OrderedPoints,ReconstPoints,FaultOrd,FaultGrad] = ...
                      FaultDetection(X,Y,fX,h,PUpar,Cpar,Hpar,SmoothPar);
%% Errors
disp('RMSE errors for ordinary fault ...')
% How much the detected points are close to the exact fault?
Gamma = @(t) [0.7+0.2*sin(2*pi*t) t];
tt = linspace(0,1,500)';
xx = Gamma(tt);
D = DistMat(xx, FaultOrd);
MinD2 = min(D.^2);
Error = sqrt(sum(MinD2)/size(FaultOrd,1));
display("RMSE distance between detected points and the exact ordinary fault is"+" "+Error);

% How much the approximate curve is close to the exact fault?
xp = ReconstPoints{1}(:,1); yp = ReconstPoints{1}(:,2);
D = DistMat(xx,[xp yp]);
MinD2 = min(D.^2);
Error = sqrt(sum(MinD2)/size([xp yp],1));
display("RMSE distance between the approximate curve and the exact ordinary fault is "+" "+Error);

disp('RMSE errors for gradient fault ...')
% How much the detected points are close to the exact fault?
Gamma = @(t) [0.4+0.1*sin(2*pi*t) t];
tt = linspace(0,1,500)';
xx = Gamma(tt);
D = DistMat(xx, FaultGrad);
MinD2 = min(D.^2);
Error = sqrt(sum(MinD2)/size(FaultGrad,1));
display("RMSE distance between detected points and the exact gradient fault is"+" "+Error);

% How much the approximate curve is close to the exact fault?
xp = ReconstPoints{2}(:,1); yp = ReconstPoints{2}(:,2);
D = DistMat(xx,[xp yp]);
MinD2 = min(D.^2);
Error = sqrt(sum(MinD2)/size([xp yp],1));
display("RMSE distance between the approximate curve and the exact gradient fault is "+" "+Error);

%% Plots
disp('Plotting ...')
figure('Name','Example 1','NumberTitle','off')  
[x,y] = meshgrid(0:0.005:1,0:0.005:1);
z = f(x,y);
surf(x,y,z);
shading interp
set(gca,'TickLabelInterpreter','latex')
set(gcf, 'Position', [300 300 400 400])
set(gca, 'XTick', [0 0.5 1])
set(gca, 'YTick', [0 0.5 1])

% detected points
figure('Name','Example 1','NumberTitle','off');
plot(XnearFault(:,1),XnearFault(:,2),'r.');
leg = legend('Detected fault points');
set(gca,'TickLabelInterpreter','latex')
set(gcf, 'Position', [300 300 400 400])
set(gca, 'XTick', [0 0.5 1]); set(gca, 'YTick', [0 0.5 1]);
set(leg,'Interpreter','latex');
xlim([0,1]); ylim([0,1.2]); box on;

% narrowed points
figure('Name','Example 1','NumberTitle','off'); hold on
p2 = plot(OrderedPoints(:,1),OrderedPoints(:,2),'ro');
p1 = plot(NarrowPoints(:,1),NarrowPoints(:,2),'b.');
set(gca,'TickLabelInterpreter','latex')
set(gcf, 'Position', [300 300 400 400])
set(gca, 'XTick', [0 0.5 1]); set(gca, 'YTick', [0 0.5 1]);
leg = legend([p1,p2],{'Narrowed points','Ordered points'});
set(leg,'Interpreter','latex');
xlim([0,1]); ylim([0,1.2]); box on;

% Plots of exact and reconstructed curves
figure('Name','Example 1','NumberTitle','off') ; hold on
gam1 = @(x,y) x-0.7-0.2*sin(2*pi*y);
gam2 = @(x,y) x-0.4-0.1*sin(2*pi*y);
p3 = fimplicit(gam1,[0 1 0 1],'k--');
fimplicit(gam2,[0 1 0 1],'k--');
for k=1:size(ReconstPoints,2)
    p2 = plot(ReconstPoints{k}(:,1),ReconstPoints{k}(:,2),'r-','LineWidth',1.5);
end
set(gca,'TickLabelInterpreter','latex')
set(gcf, 'Position', [300 300 400 400])
set(gca, 'XTick', [0 0.5 1]);set(gca, 'YTick', [0 0.5 1]);
leg = legend([p2,p3],{'Reconstructed curves','Exact curves'});
set(leg,'Interpreter','latex');
xlim([0,1]); ylim([0,1.2]); box on

%%