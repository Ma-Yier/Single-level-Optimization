function [lessMatrixLeft,lessMatrixRight,equalMatrixLeft,equalMatrixRight,objFunction] = constructMatrix(model)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明


% obtain parsing GPR info
[parInfo,nRxn,nGen,nAux,nRelation]=parseGPR(model);

% init matrix
nMet=size(model.S,1);
gprMatrixLeft=zeros(nRelation,2*nRxn+nGen+nAux);
gprMatrixRight=zeros(nRelation,1);
boundMatrixLeft=[eye(nRxn),-diag(model.ub),zeros(nRxn,nGen+nAux);-eye(nRxn),diag(model.lb),zeros(nRxn,nGen+nAux)];
boundMatrixRight=zeros(2*nRxn,1);
equalMatrixLeft=[model.S,zeros(nMet,nRxn+nGen+nAux)];
equalMatrixRight=zeros(nMet,1);
rowPoint=1;

% construct GPR matrix
for i=1:size(parInfo,1)
    
   % get each reaction GPR parsing info 
   parInfoUnit=parInfo{i,1};
   
   % has GPR
   if ~isempty(parInfoUnit{1,1})
       maxUnit=size(parInfoUnit,1);
       
       % single gene control
       if (maxUnit==1)&isempty(parInfoUnit{1,2})
           geneID=find(contains(model.genes,parInfoUnit{1,3}));
           gprMatrixLeft(rowPoint,nRxn+i)=1;
           gprMatrixLeft(rowPoint,2*nRxn+nAux+geneID)=-1;
           rowPoint=rowPoint+1;
           gprMatrixLeft(rowPoint,nRxn+i)=-1;
           gprMatrixLeft(rowPoint,2*nRxn+nAux+geneID)=1;
           rowPoint=rowPoint+1;      

       % multiple gene control
       else

           % iteratly get each auxilary GPR
           for j=1:maxUnit
               eachAuxGPR=parInfoUnit(j,:);
               logicRelation=eachAuxGPR{1,2};
               n=size(eachAuxGPR{1,3},1);
               eachAuxGPR_ctrlGene=eachAuxGPR{1,3};
               switch logicRelation
                   case 'and'

                       % obj gene
                       objLoc=locateGene(nRxn,nAux,eachAuxGPR{1,1},model);
                       gprMatrixLeft(rowPoint,objLoc)=-1;
                       gprMatrixLeft(rowPoint+1,objLoc)=n;
                       gprMatrixRight(rowPoint,1)=n-1;

                       % ctrl gene
                       for k=1:n
                           objLoc=locateGene(nRxn,nAux,eachAuxGPR_ctrlGene{k,1},model);
                           gprMatrixLeft(rowPoint,objLoc)=1;
                           gprMatrixLeft(rowPoint+1,objLoc)=-1;
                       end
                       rowPoint=rowPoint+2;

                   case 'or'

                       % obj gene
                       objLoc=locateGene(nRxn,nAux,eachAuxGPR{1,1},model);
                       gprMatrixLeft(rowPoint,objLoc)=-n;
                       gprMatrixLeft(rowPoint+1,objLoc)=1;

                       % ctrl gene
                       for k=1:n
                           objLoc=locateGene(nRxn,nAux,eachAuxGPR_ctrlGene{k,1},model);
                           gprMatrixLeft(rowPoint,objLoc)=1;
                           gprMatrixLeft(rowPoint+1,objLoc)=-1;
                       end
                       rowPoint=rowPoint+2;

                   otherwise
                       error('parsing GPR error')
               end   
           end
       end    
   end
end

% construct matrix
lessMatrixLeft=[boundMatrixLeft;gprMatrixLeft];
lessMatrixRight=[boundMatrixRight;gprMatrixRight];
objFunction=[transpose(model.c),zeros(1,nRxn+nAux+nGen)];

% end funciton
end

