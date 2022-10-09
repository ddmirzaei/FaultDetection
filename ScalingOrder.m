function ScaleOrd = ScalingOrder(op)
% This function returns the scaling order of operator op
% Inputs:
%  op: operator
% Outputs:
%  ScaleOrd: the scaling order
%
numop = length(op);
ScaleOrd = zeros(1,numop);
for k = 1:numop
    if strcmp(op{k},'1')
        ScaleOrd(k) = 0;
    elseif (strcmp(op{k},'x')||strcmp(op{k},'y')||strcmp(op{k},'z'))
        ScaleOrd(k) = 1;
    elseif (strcmp(op{k},'L')||strcmp(op{k},'xx')||strcmp(op{k},'yy')||...
            strcmp(op{k},'zz')||strcmp(op{k},'xy')||strcmp(op{k},'xz')||strcmp(op{k},'yz'))
        ScaleOrd(k) = 2;
    else   
    error('this type of op is not implemented in the ScalingOrder.m file');
    end
end
