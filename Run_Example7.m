%% Example 7 (the dolphin)
clear 
%% Function F
x = @(t) 0.001*(4/23*sin(62/33-58*t)+8/11*sin(10/9-56*t)+17/24*sin(38/35-55*t)+30/89*sin(81/23-54*t)+....
    3/17*sin(53/18-53*t)+21/38*sin(29/19-52*t)+11/35*sin(103/40-51*t)+7/16*sin(79/18-50*t)+...
    4/15*sin(270/77-49*t)+19/35*sin(59/27-48*t)+37/43*sin(71/17-47*t)+sin(18/43-45*t)+21/26*sin(37/26-44*t)+....
    27/19*sin(111/32-42*t)+8/39*sin(13/25-41*t)+23/30*sin(27/8-40*t)+23/21*sin(32/35-37*t)+18/37*sin(91/31-36*t)+...
    45/22*sin(29/37-35*t)+56/45*sin(11/8-33*t)+4/7*sin(32/19-32*t)+54/23*sin(74/29-31*t)+28/19*sin(125/33-30*t)+...
    19/9*sin(73/27-29*t)+16/17*sin(737/736-28*t)+52/33*sin(130/29-27*t)+41/23*sin(43/30-25*t)+29/20*sin(67/26-24*t)+...
    64/25*sin(136/29-23*t)+162/37*sin(59/34-21*t)+871/435*sin(199/51-20*t)+61/42*sin(58/17-19*t)+159/25*sin(77/31-17*t)+...
    241/15*sin(94/31-13*t)+259/18*sin(114/91-12*t)+356/57*sin(23/25-11*t)+2283/137*sin(23/25-10*t)+...
    1267/45*sin(139/42-9*t)+613/26*sin(41/23-8*t)+189/16*sin(122/47-6*t)+385/6*sin(151/41-5*t)+2551/38*sin(106/35-4*t)+...
    1997/18*sin(6/5-2*t)+43357/47*sin(81/26-t)-4699/35*sin(3*t+25/31)-1029/34*sin(7*t+20/21)-250/17*sin(14*t+7/40)-...
    140/17*sin(15*t+14/25)-194/29*sin(16*t+29/44)-277/52*sin(18*t+37/53)-94/41*sin(22*t+33/31)-57/28*sin(26*t+44/45)-...
    128/61*sin(34*t+11/14)-111/95*sin(38*t+55/37)-85/71*sin(39*t+4/45)-25/29*sin(43*t+129/103)-7/37*sin(46*t+9/20)...
    -17/32*sin(57*t+11/28)-5/16*sin(59*t+32/39))+1.4;
y=@(t) 0.001*(5/11*sin(163/37-59*t)+7/22*sin(19/41-58*t)+30/41*sin(1-57*t)+37/29*sin(137/57-56*t)+....
    5/7*sin(17/6-55*t)+11/39*sin(46/45-52*t)+25/28*sin(116/83-51*t)+25/34*sin(11/20-47*t)+8/27*sin(81/41-46*t)+...
    44/39*sin(78/37-45*t)+11/25*sin(107/37-44*t)+7/20*sin(7/16-41*t)+30/31*sin(19/5-40*t)+37/27*sin(148/59-39*t)+...
    44/39*sin(17/27-38*t)+13/11*sin(7/11-37*t)+28/33*sin(119/39-36*t)+27/13*sin(244/81-35*t)+13/23*sin(113/27-34*t)+...
    47/38*sin(127/32-33*t)+155/59*sin(173/45-29*t)+105/37*sin(22/43-27*t)+106/27*sin(23/37-26*t)+97/41*sin(53/29-25*t)+...
    83/45*sin(109/31-24*t)+81/31*sin(96/29-23*t)+56/37*sin(29/10-22*t)+44/13*sin(29/19-19*t)+18/5*sin(34/31-18*t)+....
    163/51*sin(75/17-17*t)+152/31*sin(61/18-16*t)+146/19*sin(47/20-15*t)+353/35*sin(55/48-14*t)+355/28*sin(102/25-12*t)+...
    1259/63*sin(71/18-11*t)+17/35*sin(125/52-10*t)+786/23*sin(23/26-6*t)+2470/41*sin(77/30-5*t)+2329/47*sin(47/21-4*t)+...
    2527/33*sin(23/14-3*t)+9931/33*sin(51/35-2*t)-11506/19*sin(t+56/67)-2081/42*sin(7*t+9/28)-537/14*sin(8*t+3/25)-...
    278/29*sin(9*t+23/33)-107/15*sin(13*t+35/26)-56/19*sin(20*t+5/9)-5/9*sin(21*t+1/34)-17/24*sin(28*t+36/23)-...
    21/11*sin(30*t+27/37)-138/83*sin(31*t+1/7)-10/17*sin(32*t+29/48)-31/63*sin(42*t+27/28)-4/27*sin(43*t+29/43)-...
    13/24*sin(48*t+5/21)-4/7*sin(49*t+29/23)-26/77*sin(50*t+29/27)-19/14*sin(53*t+61/48)+34/25*sin(54*t+37/26))+1.1;
t = 0:0.01:2*pi;
Gamma = [x(t') y(t')]/2.5;
%% Initial point set X and function values 
PointType = 'ToyExamplePoints';
X = ScatPoints2D(0,1,PointType);  
Inside = inpolygon(X(:,1),X(:,2),Gamma(:,1),Gamma(:,2));  % Scattered points in the Dolphin area (toy example).
OutF = zeros(size(Inside,1),1);                           % f = 0
OutF(Inside) = 1;                                         % f = 1
fX = OutF;                                              % Function f for toy example.
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
SmoothPar = 1;% Smoothing parameter for the cubic smoothing spline.
C_cov = 2.5;
PUpar.hcov = C_cov*h;    % spacing distance between patch centers 
PUpar.ovlp = 1.5;        % overlapping constant (uses in the narrowing step)
%%
Y = PatchCenters(PointType,0,1,PUpar.hcov); % patch centers 
%%
disp('Detection algorithm starts ...')
[NarrowPoints,XnearFault,OrderedPoints,ReconstPoints] = ...
                      FaultDetection(X,Y,fX,h,PUpar,Cpar,Hpar,SmoothPar);

%% Plots
disp('Plotting ...')
% detected points
figure('Name','Example 7','NumberTitle','off');
plot(XnearFault(:,1),XnearFault(:,2),'r.');
leg = legend('Detected fault points');
set(gca,'TickLabelInterpreter','latex')
set(gcf, 'Position', [300 300 400 400])
set(gca, 'XTick', [0 0.5 1]); set(gca, 'YTick', [0 0.5 1]);
set(leg,'Interpreter','latex');
xlim([0,1]); ylim([0,1]); box on;

% narrowed points
figure('Name','Example 7','NumberTitle','off'); hold on
p2 = plot(OrderedPoints(:,1),OrderedPoints(:,2),'ro');
p1 = plot(NarrowPoints(:,1),NarrowPoints(:,2),'b.');
set(gca,'TickLabelInterpreter','latex')
set(gcf, 'Position', [300 300 400 400])
set(gca, 'XTick', [0 0.5 1]); set(gca, 'YTick', [0 0.5 1]);
leg = legend([p1,p2],{'Narrowed points','Ordered points'});
set(leg,'Interpreter','latex');
xlim([0,1]); ylim([0,1]); box on;

% Plots of exact and reconstructed curves
figure('Name','Example 7','NumberTitle','off') ; hold on
p3 = plot(Gamma(:,1),Gamma(:,2),'k--');
for k=1:size(ReconstPoints,2)
    p2 = plot(ReconstPoints{k}(:,1),ReconstPoints{k}(:,2),'r-','LineWidth',1.5);
end
set(gca,'TickLabelInterpreter','latex')
set(gcf, 'Position', [300 300 400 400])
set(gca, 'XTick', [0 0.5 1]);set(gca, 'YTick', [0 0.5 1]);
leg = legend([p2,p3],{'Reconstructed curves','Exact curves'});
set(leg,'Interpreter','latex');
xlim([0,1]); ylim([0,1]); box on




