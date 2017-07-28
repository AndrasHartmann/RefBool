addpath ../src/

%% Get all files to booleanize
!ls ~/Work/LCSB_EPI_CBG/Data/Input/*/*EXPRESSION.txt > x.tmp
%!ls ~/Work/LCSB_EPI_CBG/Data/Differentiation_Examples/*/*EXPRESSION.txt > x.tmp
%!ls ../*EXPRESSION.txt > x.tmp


%%Initialize background
background_genes = readtable('../ReferenceDistributions/background_genes.txt', ...,
      'delimiter','\t', 'ReadVariableNames',false );

gNames = background_genes.(1);

load ../ReferenceDistributions/Thresholds_CoTFcRF.mat;

%% Do the actual booleanization
fid = fopen('x.tmp','rt');
while true
     thisline = fgetl(fid);
     if ~ischar(thisline); 
        break; 
     end  %end of file


     Expr = readtable(thisline,'delimiter', '\t');
     %align expression to dists
     indices = -ones(size(gNames));
     %Gene_ID = cell(size(gNames)); 
     for j = 1:length(indices)
        expind = find(ismember(Expr.gene_id,gNames{j}));
        if ~isempty(expind)
           indices(j)=expind;
        end
     end

     %removing non-existing indices
     mydists = dists;
     toDel = (indices==-1);

     mydists(toDel)=[];

     indices(toDel)=[];


     %vector for expression
     expression = Expr.TPM(indices);
     Gene_ID = Expr{indices,1};

     %p-values
    [ Pvals_low,Pvals_up,Pvals_inter_low,Pvals_inter_up,Pvals_low_adj,Pvals_up_adj,Pvals_inter_low_adj,Pvals_inter_up_adj ] = CalculatePValues( mydists, expression); 

    % booleanizing everything
     bool_value_all = Pvals_inter_low_adj < 0.95;

     %only lower and upper 10%
     bool_value = nan(size(Pvals_inter_low_adj));

     bool_value (Pvals_low_adj < 0.1) = 0;
     bool_value (Pvals_up_adj < 0.1) = 1;
     
     output = table(Gene_ID, expression, bool_value_all, bool_value);
     background_genes.Properties.VariableNames = {'Gene_ID', 'HGNC_symbol'};
     output = innerjoin(background_genes,output);

     %Write file
     [pathstr,name,ext] = fileparts(thisline)
     writetable(output,[pathstr filesep name '_bool' ext], 'Delimiter','\t');
end
fclose(fid);
%!rm x.tmp
