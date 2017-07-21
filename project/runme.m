%prerequisites:
% download the files 
%	https://webdav-r3lab.uni.lu/public/cbg/RefBool_ReferenceDistributions/reference/Co_TF_CRF_Biomart.txt
%	https://webdav-r3lab.uni.lu/public/cbg/RefBool_ReferenceDistributions/reference/ExpressionLibrary.txt
%e.g.:
%$wget 	https://webdav-r3lab.uni.lu/public/cbg/RefBool_ReferenceDistributions/reference/Co_TF_CRF_Biomart.txt
%$wget	https://webdav-r3lab.uni.lu/public/cbg/RefBool_ReferenceDistributions/reference/ExpressionLibrary.txt

%Adding path
addpath ../src

%Initiating parallel
%p = gcp('nocreate');
%if ~isempty(p) && p.NumWorkers ~= 6
%	delete(p)
%	parpool(6)
%end

%Reading the variables
explib = readtable('ExpressionLibrary.txt', 'delimiter', '\t');

genes =readtable('Co_TF_CRF_Biomart.txt', 'delimiter', '\t');

%filter for TFs

genes.Properties.VariableNames = {'geneID', 'geneName'};
TF_explib = innerjoin(genes, explib);

%Free some space
clear explib

A = table2array(TF_explib(:,3:end));

%tic;
Th = DetermineThresholdDistributions(A, 1000, 0.001, 'AIC', false);
%toc;
