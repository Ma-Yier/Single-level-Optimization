function [rxn,gene,geneAll] = transToRev(x,model,length)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

% pres
nGene=size(model.genes,1);
nRev=size(find(model.rev>0),1);
nRxn=size(model.rxns,1);
nAux=length-2*nRxn-nGene;

% transform nonrev negative reactions
if isfield(model,'negative')
   x(model.negative,:)=-x(model.negative,:); 
end

% transform rev reactions
geneAll=x(2*nRxn+1:2*nRxn+nAux+nGene,:);
gene=x(2*nRxn+nAux+1:2*nRxn+nAux+nGene,:);
revLogic=logical(model.rev);
x(revLogic,:)=x(revLogic,:)-x(nRxn-nRev+1:nRxn,:);
rxn=x(1:nRxn-nRev,:);

% end function
end

