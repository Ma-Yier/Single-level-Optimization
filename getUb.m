function [dualUb] = getUb(leftMatrix,rightMatrix,tag)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

% pres
[m,n]=size(leftMatrix);
objFunction=zeros(n,1);
dualUb=zeros(n,1);
objFunction(1,1)=1;

% adjust to Aineq bineq
switch tag
    case 'max'
        leftMatrix=-leftMatrix;
        rightMatrix=-rightMatrix;
    case 'min'
        % pass
        
    case others
        error('tag should be max or min')
        
end

% cal dual UB
for i=1:n
   if i>1
       objFunction(i-1,1)=0;
       objFunction(i,1)=1;
   end
   [x,fval]=cplexlp(-objFunction,leftMatrix,rightMatrix);
   dualUb(i,1)=-fval;

end

% end function
end

