function [mixLeftConMatrix,mixLeftDisMatrix,mixRightMatrix] = getDual(leftConMatrix,leftDisMatrix,rightMatrix,objFunction,tag,auxub)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
% input max <=
% input min >=

% output aineq bineq


% pres
lenDis=size(leftDisMatrix,2);
lenCon=size(leftConMatrix,2);
lenZig=size(leftDisMatrix,1)*lenDis;
lenDual=size(leftConMatrix,1);

% construct all constraints 
zigMatrix=reshape(lefDisMarix,1,[]);
dualUb=getUb(tanspose(leftConMatrix),objFunction,tag);

switch tag  % ensure min or max
    case 'max'
        
        % construct min_obj = max_obj
        mixLeftConMatrixOne=[transpose(objFunction),zigMatrix,-transpose(rightMatrix),1; ...,
                             -transpose(objFunction),-zigMatrix,transpose(rightMatrix),-1];
        mixLeftDisMatrixOne=zeros(2,lenDis);
        mixRightMatrixOne=[0;0];
        
        % primal constraints
        mixLeftDisMatrixTwo=leftDisMatrix;
        mixLeftConMatrixTwo=[leftConMatrix,zeros(size(leftConMatrix,1),lenZig+lenDual+1)];
        mixRightMatrixTwo=rightMatrix;
        
        % daul constraints
        mixLeftDisMatrixThree=zeros(size(leftConMatrix,2),lenDis);
        mixLeftConMatrixThree=[zeros(size(leftConMatrix,2),lenCon+lenZig),-transpose(leftConMatrix),zeros(size(leftConMatrix,2),1)];
        mixRightMatrixThree=-objFunction;
        
        % zig constraints
        
            % zigone 
        zigDisOne=zeros(lenZig,lenDis);
        zigConOne=[zeros(lenZig,lenCon),-eye(lenZig),zeros(lenZig,lenDual+1)];
        zigRightOne=zeros(lenZig,1);
        
            % zigtwo
        zigDisTwo=zeros(lenZig,lenDis);
        zigConTwo=[zeros(lenZig,lenCon),eye(lenZig),repmat(-eye(lenDual),lenDis,1),zeros(lenZig,1)];
        zigRightTwo=zeros(lenZig,1);
            
            % zigthree
        zigDisThree=zeros(lenZig,lenDis);      
        for i=1:lenDis
            zigDisThree((i-1)*lenDual+1:(i-1)*lenDual+lenDual,i)=-dualUb;
        end 
        zigConThree=[zeros(lenZig,lenCon),eye(lenZig),zeros(lenZig,lenDual+1)];
        zigRightThree=zeros(lenZig,1);

            % zigfour
        zigDisFour=zeros(lenZig,lenDis);
        for i=1:lenDis
            zigDisThree((i-1)*lenDual+1:(i-1)*lenDual+lenDual,i)=dualUb;
        end 
        zigConFour=[zeros(lenZig,lenCon),-eye(lenZig),repmat(eye(lenDual),lenDis,1),zeros(lenZig,1)];
        zigRightFour=repmat(dualUb,lenDis,1);
        
            % merge zig constraints
        mixLeftDisMatrixFour=[zigDisOne;zigDisTwo;zigDisThree;zigDisFour];
        mixLeftConMatrixFour=[zigConOne;zigConTwo;zigConThree;zigConFour];
        mixRightMatrixFour=[zigRightOne;zigRightTwo;zigRightThree;zigRightFour];
        
        % aux var constraints
        mixLeftDisMatrixFive=zeros(2,lenDis);
        mixLeftConMatrixFive=[zeros(2,lenCon+lenZig+lenDual),[1;-1]];
        mixRightMatrixFive=[auxub;0];

  
    case 'min'
        
        % construct min_obj = max_obj
        mixLeftConMatrixOne=[transpose(objFunction),zigMatrix,-transpose(rightMatrix),-1; ...,
                             -transpose(objFunction),-zigMatrix,transpose(rightMatrix),1];
        mixLeftDisMatrixOne=zeros(2,lenDis);
        mixRightMatrixOne=[0;0];
            
        % primal constraints
        mixLeftDisMatrixTwo=-leftDisMatrix;
        mixLeftConMatrixTwo=[-leftConMatrix,zeros(size(leftConMatrix,1),lenZig+lenDual+1)];
        mixRightMatrixTwo=-rightMatrix;
        
        % daul constraints
        mixLeftDisMatrixThree=zeros(size(leftConMatrix,2),lenDis);
        mixLeftConMatrixThree=[zeros(size(leftConMatrix,2),lenCon+lenZig),transpose(leftConMatrix),zeros(size(leftConMatrix,2),1)];
        mixRightMatrixThree=objFunction;
        
        % zig constraints
        
            % zigone 
        zigDisOne=zeros(lenZig,lenDis);
        zigConOne=[zeros(lenZig,lenCon),-eye(lenZig),zeros(lenZig,lenDual+1)];
        zigRightOne=zeros(lenZig,1);
        
            % zigtwo
        zigDisTwo=zeros(lenZig,lenDis);
        zigConTwo=[zeros(lenZig,lenCon),eye(lenZig),repmat(-eye(lenDual),lenDis,1),zeros(lenZig,1)];
        zigRightTwo=zeros(lenZig,1);
            
            % zigthree
        zigDisThree=zeros(lenZig,lenDis);      
        for i=1:lenDis
            zigDisThree((i-1)*lenDual+1:(i-1)*lenDual+lenDual,i)=-dualUb;
        end 
        zigConThree=[zeros(lenZig,lenCon),eye(lenZig),zeros(lenZig,lenDual+1)];
        zigRightThree=zeros(lenZig,1);

            % zigfour
        zigDisFour=zeros(lenZig,lenDis);
        for i=1:lenDis
            zigDisThree((i-1)*lenDual+1:(i-1)*lenDual+lenDual,i)=dualUb;
        end 
        zigConFour=[zeros(lenZig,lenCon),-eye(lenZig),repmat(eye(lenDual),lenDis,1),zeros(lenZig,1)];
        zigRightFour=repmat(dualUb,lenDis,1);
        
            % merge zig constraints
        mixLeftDisMatrixFour=[zigDisOne;zigDisTwo;zigDisThree;zigDisFour];
        mixLeftConMatrixFour=[zigConOne;zigConTwo;zigConThree;zigConFour];
        mixRightMatrixFour=[zigRightOne;zigRightTwo;zigRightThree;zigRightFour];
        
        % aux var constraints
        mixLeftDisMatrixFive=zeros(2,lenDis);
        mixLeftConMatrixFive=[zeros(2,lenCon+lenZig+lenDual),[1;-1]];
        mixRightMatrixFive=[auxub;0];
        
    case others
        error('tag should be min or max')
end

mixLeftConMatrix=[mixLeftConMatrixOne;mixLeftConMatrixTwo;mixLeftConMatrixThree;mixLeftConMatrixFour;mixLeftConMatrixFive];
mixLeftDisMatrix=[mixLeftDisMatrixOne;mixLeftDisMatrixTwo;mixLeftDisMatrixThree;mixLeftDisMatrixFour;mixLeftDisMatrixFive];
mixRightMatrix=[mixRightMatrixOne;mixRightMatrixTwo;mixRightMatrixThree;mixRightMatrixFour;mixRightMatrixFive];


% end function
end

