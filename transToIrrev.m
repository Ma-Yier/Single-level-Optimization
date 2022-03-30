function [newModel] = transToIrrev(model)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

% pres
nRxn=size(model.grRules,1);
nRevRxn=size(find(model.rev),1);
newModel=model;
logicalIndex=(newModel.rev==1);

% transfer dataset
newModel.grRules(nRxn+1:nRxn+nRevRxn,:)=model.grRules(logicalIndex,:);
newModel.rxns(nRxn+1:nRxn+nRevRxn,:)=model.rxns(logicalIndex,:);
if isfield(model,'rxnNames')
    newModel.rxnNames(nRxn+1:nRxn+nRevRxn,:)=model.rxnNames(logicalIndex,:);
end
newModel.S=[model.S, -1*model.S(:,logicalIndex)];
newModel.lb=[model.lb;zeros(nRevRxn,1)];
newModel.lb(logicalIndex,:)=0;
ubRevRxn=-1*model.lb(logicalIndex,:);
newModel.ub=[model.ub;ubRevRxn];
newModel.c=[model.c;zeros(nRevRxn,1)];
newModel.description='new_no_reversible_reaction_network';

% check if there are negative lb rxns
negativeLb=find((newModel.lb<0),1);
if ~isempty(negativeLb)
   newModel.negative=negativeLb;
   newModel.S(:,negativeLb)=-1*newModel.S(:,negativeLb);
   auxLb=newModel.ub(negativeLb,:);
   auxUb=newModel.lb(negativeLb,:);
   newModel.lb(negativeLb,:)=-1*auxLb;
   newModel.ub(negativeLb,:)=-1*auxUb;
end

% end function
end

