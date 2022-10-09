function [OrderedPoints,FaultLabel,NumFaults,FaultCycle] = OrderedSubset(X,H1,H2)
% This function chooses some sets of ordered faults points among a clode of
% points, and distinguishes between different faults
% Input:
%  X: A set of points which are narrow enough 
%  H1: Approximate distance between two consecutive ordered points
%  H2: The distance that distinguished two different faults 
% Outputs:
%  OrderedPoints: Ordered points on all faults
%  FaultLabel: label index of fault points in OrderedPoints
%  NumFaults: number of faults 
%  FaultCycle: Indices of close faults
%% Finding ordered points
OrderedPoints = [];           
OrderSet = [];                
Pos = floor(size(X,1)/2)+1;     % Index of the starting point.
FaultCount = 1;                 % a counting number  
NumFaults = 0;                  % Number of faults.
FaultCycle = [];                % Indices of close faults
while ~isempty(Pos)
    NumFaults = NumFaults+1;
    Zn = []; Zp=[];                    % Ordered points from two sides (Zn,Zp)(negative and positive sides).
    z = X(Pos,:);                      % The starting point 
    IndFz = DistMat(X,z) < H1;         % A neighborhood around point z.
    Fz = X(IndFz,:);                   % Points in a ball of radius H1 and center z
    uz = FindDirection(Fz);            % Finding the principal direction 
    Fzv = (Fz-z)*uz';                      
    [Val_zp,Pos_zp] = max(Fzv); vp = uz;  
    [Val_zn,Pos_zn] = min(Fzv); vm = uz;  
    znew = Fz(Pos_zn(1),:);            
    FzSave = Fz;                       
    cycle = 0;                         
    OrderSet = [OrderSet;z];         
    if FaultCount == 1, DistFtoZnew = H1; else, [DistFtoZnew,IndDist] = min(DistMat(znew,OrderSet));end 
    if (FaultCount~=1) && isequal(OrderSet(IndDist,:),z)
       DistFtoZnew = H1;              
    end 
    %  Finding Zn part of the fault
    while sign(Val_zn(1)) < 0 && DistFtoZnew(1) > H1/2
        Zn = [znew;Zn];                
        OrderSet = [OrderSet;znew];   
        IndFz = DistMat(X,znew)<H1;   
        Fz = X(IndFz,:);              
        uz = FindDirection(Fz);       
        if dot(vm,uz) < 0, uz = -uz; end 
        vm = uz;
        Fzv = (Fz-znew)*uz';
        [Val_zn,Pos_zn] = min(Fzv);    
        znew = Fz(Pos_zn(1),:);        
        [DistFtoZnew,IndDist] = min(DistMat(znew,OrderSet(1:end-1,:)));
    end
    % Check if part Zn of the fault is a cycle (intersects itself).
    if sign(Val_zn(1)) < 0 && DistFtoZnew(1) <= H1/2  
        if ismember(OrderSet(IndDist,:),[z;Zn])       
             FaultCycle = [FaultCycle,NumFaults];
             OS1 = OrderSet(DistMat(OrderSet,znew)<H1,:);
             OS2 = OS1(DistMat(OS1,OrderSet(end,:))>0,:);
             [val,BestIndz] = min(DistMat(OS2,OrderSet(end,:)));
             Zn = [OS2(BestIndz,:);Zn];                          
             if distance(znew,z)<H1/2, cycle = 1; end            
         end
    end
    %  Finding Zp part of the fault
    Fz = FzSave;                       
    znew = Fz(Pos_zp(1),:);            
    [DistFtoZnew,IndDist] = min(DistMat(znew,OrderSet));
    if OrderSet(IndDist,:) == z, DistFtoZnew = H1; end
    while sign(Val_zp(1)) > 0 && DistFtoZnew(1) > H1/2
        if cycle == 1, break; end          
        Zp = [Zp;znew];                   
        OrderSet = [OrderSet;znew];     
        IndFz = DistMat(X,znew)<H1;     
        Fz = X(IndFz,:);                
        uz = FindDirection(Fz);         
        if dot(vp,uz) < 0, uz = -uz; end 
        vp = uz;
        Fzv = (Fz-znew)*uz';
        [Val_zp,Pos_zp] = max(Fzv);      
        znew = Fz(Pos_zp(1),:);          
        [DistFtoZnew,IndDist] = min(DistMat(znew,OrderSet(1:end-1,:)));
    end
    % Check if Zp part intersects [Zn;z;Zp] and creates a cycle
    if sign(Val_zp(1)) > 0 && DistFtoZnew(1) <= H1/2
        if ismember(OrderSet(IndDist,:),[Zn;z;Zp])
            if ~isempty(Zp)                
                OS1 = OrderSet(DistMat(OrderSet,znew)<H1,:);
                OS2 = OS1(DistMat(OS1,OrderSet(end,:))>0,:);
                [val,BestIndz] = min(DistMat(OS2,OrderSet(end,:)));
                Zp=[Zp;OS2(BestIndz,:)];
            elseif isempty(Zp)&& cycle~=1  
                OS1 = OrderSet(DistMat(OrderSet,znew)<H1,:);
                OS2 = OS1(DistMat(OS1,z)>0,:);
                [val,BestIndz] = min(DistMat(OS2,z));
                Zp = [Zp;OS2(BestIndz,:)];
             end
        end
    end
    % Check if there exists another fault 
    OrderedPoints = [OrderedPoints;Zn;z;Zp];     
    FaultLabel(FaultCount) = size(OrderedPoints,1);
    DistAll = DistMat(OrderedPoints,X);    
    minval = min(DistAll);
    IndRem = find(minval > H2);  
    if ~isempty(IndRem), Pos = IndRem(1); else, Pos = []; end
    FaultCount = FaultCount+1;
end
%%
