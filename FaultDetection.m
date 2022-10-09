function [NarrowPoints,XnearFault,OrderedPoints,ReconstPoints,FaultOrd,FaultGrad]=...
              FaultDetection(X,Y,fX,h,PUpar,Cpar,Hpar,SmoothPar)
% This function (1) detects a set of fault points around the fault curves
%               (2) narrows the the detected points
%               (3) extracts a set of ordered points out of narrowed points 
%               (3) fixes the intersections
%               (4) classifies the fault points 
%
% Inputs:
%  X: initial points
%  Y: patch centers 
%  fX: function values
%  h: Fill distance.
%  PUpar: hcov, C_ovlp parameters
%  Cpar: CL, CG, CM, eta parameters
%  Hpar: H1,H2,H3,H4, H5 parameters
%  SmoothPar: Smoothing parameter for the cubic smoothing spline.

% Outputs:
% NarrowPoints: narrowed points 
% XnearFault: detected points near faults
% FaultOrd: points near Ordinary faults
% FaultGrad: points near Gradient faults
% OrderedPoints: ordere fault points
% ReconstPoints: {1*number of fault} cells containing points on the cubic smoothing spline recostructions
%%
global RBFtype RBFpar PolyOrder  

[IndXY,dist] = knnsearch(X,Y,'k',12);
rho = max(dist,[],2);
IndYY = PointsInPatch(Y,Y,1.5*max(rho));   

RBFtype = 'p'; RBFpar = 3;                           % Polyharmonic RBF: r^RBFpar
PolyOrder = floor(RBFpar/2)+1;                       % Polynomial order
LapApp = RBF_PU(X,Y,IndXY,IndYY,fX,rho,'L');         % RBF-PU interpolant for Laplacian
LapNorm = abs(LapApp);

RBFtype = 'p'; RBFpar = 1;                           % Polyharmonic spline RBF: r^RBFpar
PolyOrder = floor(RBFpar/2)+1;                       % Polynomial order
GradApp = RBF_PU(X,Y,IndXY,IndYY,fX,rho,{'x','y'});  % RBF-PU interpolant for gradient
gx = GradApp{1}; gy = GradApp{2};
GradNorm = sqrt(gx.^2+gy.^2); 

%% Laplace indicator for fault detection
CL = Cpar.L; CG = Cpar.G; CM = Cpar.M; 
alpha2 = CL/(h^1);
IdxLapAlfa2 = (LapNorm > alpha2);
LapF = LapNorm(IdxLapAlfa2); 
med = median(LapF);
delta2 = CM*med;
IdxLapDel2 = (LapNorm > delta2);
FaultOrdGrad = X(IdxLapDel2,:);                               % points near both gradient and ordinary faults

%% Gradient indicator for fault detection.
alpha1 = CG/(h^(1/2));
GradFault0 = GradNorm(IdxLapDel2);
IdxGradAlfa1 = (GradFault0 > alpha1);
GradAlfa1 = GradFault0(IdxGradAlfa1); 
med = median(GradAlfa1);
delta1 = CM*med;
IdxGradDel1 = (GradFault0 > delta1);
FaultOrd = FaultOrdGrad(IdxGradDel1,:);                       % points near ordinary faults
IdxFaultOrdGrad=find(boolean(ones(size(FaultOrdGrad,1),1)));  % Indices of points near both ordinary and gradient points
IdxFaultOrd = find(IdxGradDel1);                              % Indices of of points near ordinary fault
IdxFaultGrad = find(~IdxGradDel1);                            % Indices of points near gradient fault (corrected)

%% Deleting gradient fault points which are near ordinary fault points.
Ind = [];
for i = 1:size(IdxFaultGrad,1)
    r = find(DistMat(FaultOrdGrad(IdxFaultOrdGrad,:),FaultOrdGrad(IdxFaultGrad(i),:)) < Cpar.eta);
    CardB = length(r);
    CardBF = length(intersect(IdxFaultOrdGrad(r),IdxFaultOrd));
    if CardBF >= (1/5)*CardB
        Ind = [Ind,i];
    end
end
IdxFaultGrad(Ind) = [];
FaultGrad = FaultOrdGrad(IdxFaultGrad,:);           % Final points near gradient fault 
XnearFault = [FaultOrd;FaultGrad];          % Final points near both ordianry and gradient faults   

%% Narrowing step 
C_ovlp = PUpar.ovlp; hcov = PUpar.hcov;
PatchRad = C_ovlp*hcov;   % patch radius for PU method in the narrowing step  
NarrowPoints = XnearFault;
for n =1:2                             
    PatchCentersOnFault = OrderedSubset(NarrowPoints,PatchRad,PatchRad);
    NarrowPoints = PCR_PU(NarrowPoints,PatchCentersOnFault,PatchRad);  % PU and PCR to narrow detected points
end

%% Finding a set of ordered points out of narrowing faults points
disp('Constructing ordered points ...')
[OrderedPoints,FaultLabel,NumFaults,FaultCycle] = OrderedSubset(NarrowPoints,Hpar.H1,Hpar.H2);

%% Fixing intersections 
[OrderedPoints,FaultLabel] = FixingIntersections(OrderedPoints,FaultLabel,Hpar,NumFaults,FaultCycle);

%% Finding reconstruction curves by cubic smoothing spline.
FaultLabel = [0 FaultLabel];                                               
ReconstPoints = {};
for k=1:length(FaultLabel)-1                                        
    Fk = OrderedPoints(FaultLabel(k)+1:FaultLabel(k+1),:);                
    ReconstPoints{k} = ParCubicSpline(Fk,500,SmoothPar);      
end
