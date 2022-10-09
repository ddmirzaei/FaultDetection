function [OrdPts,Labels] = FixingIntersections(OrdPts,Labels,Hpar,NumFaults,Cycle)
% This function finds intersection points and modifies the ordered points
%
% Inputs:
%   OrdPts: Old sets of ordered points 
%   Labels: Old lablel for each fault from OrdPts
%   Hpar: Parameters H3, H4, H5  
%   NumFaults: Number of initial faults
%   Cycle: Indices of cycle faults
% Outputs:
%   OrdPts: New sets of ordered points 
%   Labels: New Lablel for each fault from OrdPts
%% Creating an empty information cell for the intersection points of diffrerent faults.
H3 = Hpar.H3; H4 = Hpar.H4; H5 = Hpar.H5; 
Intersect = cell(1,NumFaults);   
for i = 1:NumFaults              
    Intersect{i} = cell(1,3);      
end                              
%% (steps 1, 2 and 3): 
% Extending head (first two points of the fault) and tail (two end points  of the fault) for each non cycle fault
IndFaults = setdiff(1:NumFaults,Cycle);  
Labels=[0,Labels];                                    
for k = IndFaults    
    % steps 1 and 2
    ExtHead = OrdPts(Labels(k)+1,:)+(H3/norm(OrdPts(Labels(k)+1,:)-...
              OrdPts(Labels(k)+2,:),2))*(OrdPts(Labels(k)+1,:)-OrdPts(Labels(k)+2,:));
    ExtTail = OrdPts(Labels(k+1),:)+(H3/norm(OrdPts(Labels(k+1),:)-...
              OrdPts(Labels(k+1)-1,:),2))*(OrdPts(Labels(k+1),:)-OrdPts(Labels(k+1)-1,:));
    LineH = [OrdPts(Labels(k)+1,:);ExtHead];                     % Adding end point to the head 
    LineT = [OrdPts(Labels(k+1),:);ExtTail];                     % Adding end point to the tail 
    RemPts = OrdPts([1:Labels(k),Labels(k+1)+1:end],:);          % Deleting the kth fault from ordered points
    % Head and tail formulas 
    ah = (LineH(2,2)-LineH(1,2))/(LineH(2,1)-LineH(1,1));  
    bh = LineH(1,2)-ah*LineH(1,1);                         
    at = (LineT(2,2)-LineT(1,2))/(LineT(2,1)-LineT(1,1));  
    bt = LineT(1,2)-at*LineT(1,1);                         
    % step 3 
    IndxP = ((ah*RemPts(:,1)+bh) >= RemPts(:,2));
    PtsOnPosH = RemPts(IndxP,:);                % Ordered points on the positive side of the head
    PtsOnNegH = RemPts(~IndxP,:);               % Ordered points on the negative side of the head
    IndxN=((at*RemPts(:,1)+bt) >= RemPts(:,2));
    PtsOnPosT = RemPts(IndxN,:);                % Ordered points on the positive side of the tail
    PtsOnNegT = RemPts(~IndxN,:);               % Ordered points on the negative side of the tail
    % Finding Closest ordered points on the positive and negative side of the the head and tail 
    zpH = []; znH = []; zpT = []; znT = [];              
    if ~isempty(PtsOnPosH)                               
      [zpHdist,zpHind] = min(DistMat(PtsOnPosH,ExtHead));
      if zpHdist <= H4, zpH = PtsOnPosH(zpHind(1),:); end 
    end
    if ~isempty(PtsOnNegH)
      [znHdist,znHind] = min(DistMat(PtsOnNegH,ExtHead));
      if znHdist <= H4, znH = PtsOnNegH(znHind(1),:); end 
    end
    if ~isempty(PtsOnPosT)
      [zpTdist,zpTind] = min(DistMat(PtsOnPosT,ExtTail));
      if zpTdist <= H4, zpT = PtsOnPosT(zpTind(1),:); end 
    end
    if ~isempty(PtsOnNegT)
      [znTdist,znTind] = min(DistMat(PtsOnNegT,ExtTail));
      if znTdist <= H4, znT = PtsOnNegT(znTind(1),:); end 
    end
    % Determining the intersection point of the extension of tail or head
    %    with the line segment between zp and zn
    if (~isempty(zpH) && ~isempty(znH))        
        a = (znH(2)-zpH(2))/(znH(1)-zpH(1));   
        b = zpH(2)-a*zpH(1);                
        x = (b-bh)/(ah-a);      
        y = (ah*b-bh*a)/(ah-a);
        Intersect{k}{1}=k; Intersect{k}{2} = [x,y];
    elseif ~isempty(zpH)                              
        Intersect{k}{1} = k; Intersect{k}{2} = zpH;   
    elseif ~isempty(znH)                              
        Intersect{k}{1} = k; Intersect{k}{2} = znH;   
    end
    if (~isempty(zpT) && ~isempty(znT))               
        a = (znT(2)-zpT(2))/(znT(1)-zpT(1));
        b = zpT(2)-a*zpT(1);
        x = (b-bt)/(at-a);     
        y = (at*b-bt*a)/(at-a);
        Intersect{k}{1} = k; Intersect{k}{3} = [x,y];
    elseif ~isempty(zpT)                              
        Intersect{k}{1} = k; Intersect{k}{3} = zpT;       
    elseif ~isempty(znT)                              
        Intersect{k}{1} = k; Intersect{k}{3} = znT;       
    end
end
% whether the intersection point is derived from head or tail of a fault
HeadTail = [];                             
for i = 1:NumFaults                         
    if ~isempty(Intersect{i}{1})          
        if ~isempty(Intersect{i}{2})     
            HeadTail=[HeadTail;i,2];    
        end
        if ~isempty(Intersect{i}{3})
            HeadTail=[HeadTail;i,3];    
        end
    end
end

% step 4: finding intersection points whose their distances are lese than H5
count = 1;      
NearbyIntersect = {};
while ~isempty(HeadTail)
    row = [];
    NearInd = [];
    for i=1:size(HeadTail,1) 
        if DistMat(Intersect{HeadTail(1,1)}{HeadTail(1,2)},Intersect{HeadTail(i,1)}{HeadTail(i,2)}) < H5
            NearInd = [NearInd;HeadTail(i,:)];
            row = [row,i];
        end
    end
    NearbyIntersect{count} = NearInd;      
    HeadTail(row,:) = [];
    count = count+1;
end
% close points are replaced by their average
AvgIntersect=[];                        
for i = 1:length(NearbyIntersect)          
    NearPts = [];
    for j = 1:size(NearbyIntersect{i},1)
        NearPts = [NearPts;Intersect{NearbyIntersect{i}(j,1)}{NearbyIntersect{i}(j,2)}];
    end
    AvgIntersect(i,:) = mean(NearPts,1);
end

% Step 5: constructing the new oredered sets by adding intersections
NewOrd = {}; 
for k = 1:length(Labels)-1
    NewOrd{k} = OrdPts(Labels(k)+1:Labels(k+1),:);
end
for i = 1:length(NearbyIntersect)                     
    if size(NearbyIntersect{i},1) == 1 
        IndxFault = NearbyIntersect{i}(1,1);
        if NearbyIntersect{i}(1,2) == 2  
            NewOrd{IndxFault} = [AvgIntersect(i,:);NewOrd{IndxFault}];
        else                        
            NewOrd{IndxFault} = [NewOrd{IndxFault};AvgIntersect(i,:)]; 
        end
    elseif size(NearbyIntersect{i},1) > 1    
        KindIntrs1 = NearbyIntersect{i}(1,2);
        KindIntrs2 = NearbyIntersect{i}(2,2);
        if ((KindIntrs1 == 2) && (KindIntrs2 == 2))
            FirstFault = min(NearbyIntersect{i}(1,1),NearbyIntersect{i}(2,1));    
            SecondFault = max(NearbyIntersect{i}(1,1),NearbyIntersect{i}(2,1));  
            NewOrd{FirstFault} = [flip(NewOrd{SecondFault},1);AvgIntersect(i,:);NewOrd{FirstFault}]; 
            NewOrd{SecondFault} = [];                 
        elseif ((KindIntrs1 == 3) && (KindIntrs2 == 3)) 
            FirstFault = min(NearbyIntersect{i}(1,1),NearbyIntersect{i}(2,1));    
            SecondFault = max(NearbyIntersect{i}(1,1),NearbyIntersect{i}(2,1));   
            NewOrd{FirstFault} = [NewOrd{FirstFault};AvgIntersect(i,:);flip(NewOrd{SecondFault},1)];
            NewOrd{SecondFault}=[];
        else             
            [t,IndH] = max([NearbyIntersect{i}(1,2),NearbyIntersect{i}(2,2)]);
            [h,IndT] = min([NearbyIntersect{i}(1,2),NearbyIntersect{i}(2,2)]);
            FirstFault = NearbyIntersect{i}(IndH,1);
            SecondFault = NearbyIntersect{i}(IndT,1);                     
            NewOrd{FirstFault} = [NewOrd{FirstFault};AvgIntersect(i,:);NewOrd{SecondFault}];
            NewOrd{SecondFault} = [];
        end
    end
end
OrdPts = []; Labels = []; num = 0; count = 1;
for i = 1:length(NewOrd)
    if ~isempty(NewOrd{i})
        Labels(count) = [num+size(NewOrd{i},1)];
        num = Labels(count);
        OrdPts = [OrdPts;NewOrd{i}];
        count = count+1;
    end
end
%%

