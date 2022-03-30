function [objLoc] = locateGene(nRxn,nAux,objGene,model)
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

% get name length
geneNameLen=size(objGene,2);

% if name is aux.. real.. or 
if geneNameLen>2
    
    % get first 3 alpha
    geneNameHead=objGene(1:3);
    switch geneNameHead

        % real case 
        case 'rea'
            objLoc=str2num(objGene(5))+nRxn;

        % aux case
        case 'aux'
            objLoc=str2num(objGene(4))+2*nRxn;

        % gene case
        otherwise
            objLoc=find(contains(model.genes,objGene))+2*nRxn+nAux;

    end
    
else
    
    % gene case
    objLoc=find(contains(model.genes,objGene))+2*nRxn+nAux;    
end

% end fucntion
end

