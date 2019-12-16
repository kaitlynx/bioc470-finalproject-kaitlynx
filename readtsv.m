function [readcount,cluster,genes] = readtsv(filename)
%READTSV Reads a tsv file and returns the read count matrix and cluster
%labeling
opts = detectImportOptions(filename,'FileType','text');
opts.VariableTypes(2) = {'char'};
tab = readtable(filename,opts);
readcount = table2array(tab(:,3:size(tab,2)));
cluster = table2array(tab(:,2));
genes = tab(:,3:size(tab,2)).Properties.VariableNames;
end
