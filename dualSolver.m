function [knockOut,gene,fval,sr] = dualSolver(model,biomassID,targetID,inputID,inputMass,alpha,beta)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明

% alpha and beta should be 0~1
if (alpha>1)|(beta>1)
   error('alpha or beta value error');   
end

% init parameters
model.lb(inputID)=inputMass;
model.ub(inputID)=inputMass;

% obtain TMPR and TMGR
xTMGR=cplexlp(-model.c,[],[],model.S,model.b,model.lb,model.ub);
TMGR=xTMGR(biomassID,1);
model.c(biomassID,1)=0;
model.c(targetID,1)=1;
xTMPR=cplexlp(-model.c,[],[],model.S,model.b,model.lb,model.ub);
TMPR=xTMPR(targetID,1);
model.lb(biomassID)=0.1*TMGR;

% prepare for dualSolver
newModel=transToIrrev(model);
[lessMatrixLeft,lessMatrixRight,equalMatrixLeft,equalMatrixRight,objFunction] ...,
                                                   =constructMatrix(newModel); %max biomass
mInPrimalMatrix=2*size(equalMatrixLeft,1)+size(lessMatrixLeft,1);
nInPrimalMatrix=size(equalMatrixLeft,2);

% normalize (max)
inPrimalLessMatrixLeft=[equalMatrixLeft;-equalMatrixLeft;lessMatrixLeft];
inPrimalLessMatrixRight=[equalMatrixRight;-equalMatrixRight;lessMatrixRight];
inPrimalObjFunction=objFunction;

% get inner dual (min)
inDualGreatMatrixLeft=transpose(inPrimalLessMatrixLeft);
inDualGreatMatrixRight=transpose(objFunction);
inDualObjFunction=transpose(inPrimalLessMatrixRight);
[mInDualMatrix,nInDualMatrix]=size(inDualGreatMatrixLeft);

% eliminate inner max/min and construct outer optimization (min)
outPrimalObjFunction=zeros(1,nInPrimalMatrix+nInDualMatrix+1);
outPrimalObjFunction(1,targetID)=1;
outPrimalGreatMatrixRight=[zeros(2,1); ...,
                           -inPrimalLessMatrixRight; ...,
                           inDualGreatMatrixRight; ...,
                           -1*alpha*TMGR; ...,
                           0];
outPrimalGreatMatrixLeft=[-inPrimalObjFunction,inDualObjFunction,-1; ...,
                          inPrimalObjFunction,-inDualObjFunction,1; ...,
                          -inPrimalLessMatrixLeft,zeros(mInPrimalMatrix,nInDualMatrix+1); ...,
                          zeros(mInDualMatrix,nInPrimalMatrix),inDualGreatMatrixLeft,zeros(mInDualMatrix,1); ...,
                          zeros(2,nInPrimalMatrix+nInDualMatrix),[-1;1]];
[mOutPrimalMatrix,nOutPrimalMatrix]=size(outPrimalGreatMatrixLeft);
                      
% get out dual (max)
outDualObjFunction=transpose(outPrimalGreatMatrixRight);
outDualLessMatrixLeft=transpose(outPrimalGreatMatrixLeft);
outDualLessMatrixRight=transpose(outPrimalObjFunction);
[mOutDualMatrix,nOutDualMatrix]=size(outDualLessMatrixLeft);

% eliminate outer min/max
singleOptMatrixLeft=[[outPrimalObjFunction;-outPrimalObjFunction],[-outDualObjFunction;outDualObjFunction],[-1;1]; ...,
                     -outPrimalGreatMatrixLeft,zeros(mOutPrimalMatrix,nOutDualMatrix+1); ...,
                     zeros(mOutDualMatrix,nOutPrimalMatrix),outDualLessMatrixLeft,zeros(mOutDualMatrix,1); ...,
                     zeros(2,nOutPrimalMatrix+nOutDualMatrix),[1;-1]];
singleOptMatrixRight=[zeros(2,1); ...,
                      -outPrimalGreatMatrixRight; ...,
                      outDualLessMatrixRight; ...,
                      beta*TMPR; ...,
                      0];
singleOptObjFunction=zeros(nOutPrimalMatrix+nOutDualMatrix+1,1);
singleOptObjFunction(targetID,1)=1;

% variable type
% replem()
ctype="";
for i=1:size(model.S,2)
   ctype=ctype+"C"; 
end
for i=1:nInPrimalMatrix-size(model.S,2)
   ctype=ctype+"B"; 
end
for i=1:(nInDualMatrix+1+nOutDualMatrix+1)
    ctype=ctype+"C"; 
end
ctype=char(ctype);

% solve
[x,fval] = cplexmilp(-singleOptObjFunction,singleOptMatrixLeft,singleOptMatrixRight,[],[],[],[],[],[],[],ctype);

% process result --> input format
[rxn,gene,geneAll]=transToRev(x,newModel,nInPrimalMatrix);

% get reaction knockout
rxnKnock = parseRxnKnock(model,geneAll);

% verify the strategy
knockOut=(rxnKnock>0);
model.lb=model.lb.*knockOut;
model.ub=model.ub.*knockOut;
model.c(biomassID,1)=1;
model.c(targetID,1)=0;
[xVerify,xBiomass]=cplexlp(-model.c,[],[],model.S,model.b,model.lb,model.ub);
xBiomass=-1*xBiomass;

% success or not
if (xBiomass>=0.1*TMGR)&(xVerify(targetID,1)>=0.1*TMPR)
    sr='success';
else
    sr='false';   
end

% end function
end

