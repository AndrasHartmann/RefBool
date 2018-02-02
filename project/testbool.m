examples = readtable('~/Work/LCSB_EPI_CBG/Data/HighLowEffi/Trans_Diff_TF_Examples_HighLow.txt', ...,
      'delimiter','\t', 'ReadVariableNames',true);

myfiles = table2array (readtable('x.txt', 'ReadVariableNames',false));

for i = 1:size(examples,1)

	upregulated = strsplit(examples.Genes_Up{i},',')

	matches = strfind(myfiles, examples.Source(i));
	source_file = myfiles{find(~cellfun(@isempty,matches))};
	source_file = regexprep(source_file,'Expression','Expression_bool');
	source = readtable(source_file, 'delimiter','\t', 'ReadVariableNames',true);

	matches = strfind(myfiles, examples.Target(i));
	target_file = myfiles{find(~cellfun(@isempty,matches))};
	target_file = regexprep(target_file,'Expression','Expression_bool');
	target = readtable(target_file, 'delimiter','\t', 'ReadVariableNames',true);

	for j = 1:length(upregulated)
		any(table2array(source(strcmp(source.HGNC_symbol,upregulated{j}),'bool_value_all')))
		any(table2array(target(strcmp(target.HGNC_symbol,upregulated{j}),'bool_value_all')))
	end
end

