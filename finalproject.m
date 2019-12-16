%% Read tsv data. There are 3 foveal samples and 3 peripheral samples. 
[fovea1_mat, fovea1_clusters,genes] = readtsv('GSM3745992_fovea_donor_1_expression.tsv');
[fovea2_mat, fovea2_clusters,~] = readtsv('GSM3745993_fovea_donor_2_expression.tsv');
[fovea3_mat, fovea3_clusters,~] = readtsv('GSM3745994_fovea_donor_3_expression.tsv');
[peripheral1_mat, peripheral1_clusters,~] = readtsv('GSM3745995_peripheral_donor_1_expression.tsv');
[peripheral2_mat, peripheral2_clusters,~] = readtsv('GSM3745996_peripheral_donor_2_expression.tsv');
[peripheral3_mat, peripheral3_clusters,~] = readtsv('GSM3745997_peripheral_donor_3_expression.tsv');
save('rawdatafovea.mat', 'fovea1_mat', 'fovea1_clusters', 'fovea2_mat', 'fovea2_clusters', 'fovea3_mat', 'fovea3_clusters','genes')
save('rawdataperipheral.mat', 'peripheral1_mat', 'peripheral1_clusters', 'peripheral2_mat', 'peripheral2_clusters', 'peripheral3_mat', 'peripheral3_clusters','genes')

%% Combine data into single matrix
alldata = [fovea1_mat;fovea2_mat;fovea3_mat;peripheral1_mat;peripheral2_mat;peripheral3_mat];
samples = ["fovea1"; 
    "fovea2"; 
    "fovea3"; 
    "peripheral1"; 
    "peripheral2"; 
    "peripheral3"];
sample = repelem(samples, [size(fovea1_mat,1);
    size(fovea2_mat,1);
    size(fovea3_mat,1);
    size(peripheral1_mat,1);
    size(peripheral2_mat,1);
    size(peripheral3_mat,1)]);
regions = ["fovea"; "peripheral"];
region = repelem(regions, [size(fovea1_mat,1) + size(fovea2_mat,1) + size(fovea3_mat,1);
    size(peripheral1_mat,1) + size(peripheral2_mat,1) + size(peripheral3_mat,1)]);
cluster = [fovea1_clusters;
    fovea2_clusters;
    fovea3_clusters;
    peripheral1_clusters;
    peripheral2_clusters;
    peripheral3_clusters];
save('combineddata.mat', 'alldata', 'sample','region','cluster');

%% Filter genes, all genes expressed in >= 3 cells
filtereddata = alldata;
filtereddata(:, ~(sum(alldata > 0) >= 3)) = [];
filteredgenes = genes;
filteredgenes(~(sum(alldata > 0) >= 3)) = [];
save('filterreddata.mat','filtereddata','filteredgenes');

%% Identify highly variable genes
x = exp(filtereddata) - 1;
mu = log(mean(x) + 1);
dispersion = log(var(x)./mean(x));
figure; plot(mu,dispersion,'r.');

% zscore
[~, ~, bin] = histcounts(mu, 20);
zsc = zeros(1, size(x,2));
for ii = 1:20
    indices = bin == ii;
    k = find(indices);
    scores = zscore(dispersion(indices));
    for jj = 1:length(scores)
        zsc(k(jj)) = scores(jj);
    end
end

% plot
figure; plot(mu,zsc,'.');
x = zsc > 0.75 & mu > 0.0125 & mu < 3;
variablegenes = filteredgenes(x);
text(mu(x), zsc(x), variablegenes,'VerticalAlignment','top','HorizontalAlignment','left','FontSize',6);
xlabel('Average Expression'); ylabel('Dispersion');

% threshold
x = zsc > 0.75 & mu > 0.0125 & mu < 3;
sum(x);
variablegenes = filteredgenes(x);
find(strcmp(variablegenes,'PDE6H'));

save('hvg.mat', 'mu','dispersion','zsc','variablegenes','x');

%% Scale data
normdata = normalize(filtereddata);

%% PCA with variable genes
rng('default');
[coeff, sc, eig] = pca(normdata(:,x));

%% Find significant PCs
plot(1:30,std(sc(:,1:30)),'.','MarkerSize',12);
xlabel('PC'); ylabel('Standard Deviation of PC');

%% t-SNE with significant PCs
rng('default');
X = tsne(sc(:,1:16));
figure;
gscatter(X(:,1),X(:,2),sample);
legend('Location','eastoutside')
xlabel('tSNE 1'); ylabel('tSNE 2');
title('tSNE plot grouped by sample');
figure;
gscatter(X(:,1),X(:,2),region);
legend('Location','eastoutside')
xlabel('tSNE 1'); ylabel('tSNE 2');
title('tSNE plot grouped by region');



%% k means clustering
rng('default');
idx = kmeans(sc(:,1:16),18);

figure;
gscatter(X(:,1),X(:,2),idx);
legend('Location','eastoutside')
xlabel('tSNE 1'); ylabel('tSNE 2');
title('tSNE plot grouped by cluster');

% Label cluster number
cluster = string(idx);
uniqueclus = string(unique(cluster));
labelxy = zeros(length(uniqueclus),2);
for ii = 1:length(uniqueclus)
    labelxy(ii,:) = mean(X(idx == str2double(uniqueclus{ii}),:));
end
text(labelxy(:,1),labelxy(:,2),uniqueclus);

save('clustering.mat', 'sc', 'X', 'idx','cluster','uniqueclus','normdata');

%% Plot PCs
hold on;
plot(sc(:,1),sc(:,2),'r.');
plot(sc(end,1),sc(end,2),'k.');
xlabel('PC1'); ylabel('PC2');
hold off;

%% Plot heatmap overlaying clusters for known cell markers
genelist = ["PDE6H" "ARR3" "RHO" "NRL" "GLUL" "RLBP1" "TRPM1" "GNG13" "SNCG" "NEFM" "FLT1" "VWF"];
figure;
for ii = 1:length(genelist)
    subplot(3,4,ii); scatter(X(:,1),X(:,2),1,alldata(:,find(strcmp(genes,genelist(ii)))));
    colormap('jet');
    title(genelist(ii));
    xlabel('tSNE 1'); ylabel('tSNE 2');
end

%% Plot gene expression by cluster
cluster = idx;
genelist = ["PDE6H" "ARR3" "RHO" "NRL" "GLUL" "RLBP1" "TRPM1" "GNG13" "SNCG" "NEFM" "FLT1" "VWF"];
figure;
for ii = 1:length(genelist)
    subplot(3,4,ii); plot(cluster,alldata(:,find(strcmp(genes,genelist(ii)))),'.');
    title(genelist(ii));
end

%% Plot boxplots by cluster
cluster = idx;
genelist = ["PDE6H" "ARR3" "RHO" "NRL" "GLUL" "RLBP1" "TRPM1" "GNG13" "SNCG" "NEFM" "FLT1" "VWF"];
figure;
for ii = 1:length(genelist)
    subplot(3,4,ii); boxplot(alldata(:,find(strcmp(genes,genelist(ii)))),cluster);
    title(genelist(ii));
end

%% Differential gene analysis
cluster = string(idx);
clustereddata = zeros(length(uniqueclus),size(filtereddata,2));
for ii = 1:length(uniqueclus)
    clustereddata(ii,:) = mean(filtereddata(strcmp(cluster,uniqueclus(ii)),:));
end

% Calculate fold change for one cluster compared to all other cells
markers = cell(length(uniqueclus)*10, 2);
for ii = 1:length(uniqueclus)
    otherclusters = filtereddata(~strcmp(cluster,uniqueclus(ii)),:);
    thiscluster = filtereddata(strcmp(cluster,uniqueclus(ii)),:);
    thisclusterfilt = thiscluster;
    otherclusterfilt = otherclusters;
    thisclusterfilt(:, ~(sum(thiscluster > 0) >= 0.25*size(thiscluster,1) | sum(otherclusters > 0) >= 0.25*size(otherclusters,1))) = [];
    otherclusterfilt(:, ~(sum(thiscluster > 0) >= 0.25*size(thiscluster,1) | sum(otherclusters > 0) >= 0.25*size(otherclusters,1))) = [];
    tempfilter = filteredgenes;
    tempfilter(:, ~(sum(thiscluster > 0) >= 0.25*size(thiscluster,1) | sum(otherclusters > 0) >= 0.25*size(otherclusters,1))) = [];
    
    foldchange = mean(exp(thisclusterfilt)) ./ mean(exp(otherclusterfilt));
    log2fc = log2(foldchange);
    [highfc, sortidx] = sort(log2fc, 'descend');
    sortedgenes = tempfilter(sortidx);
    markers((ii-1)*10+1:(ii-1)*10+10,1) = sortedgenes(1:10)';
    markers((ii-1)*10+1:(ii-1)*10+10,2) = {uniqueclus(ii)};
    disp(uniqueclus(ii));
    disp(sortedgenes(1:30));
end

% code for p-values commented out
%     pval = zeros(1,size(thisclusterfilt,2));
%     for jj = 1:size(thisclusterfilt,2)
%         [h,p] = ttest2(exp(thisclusterfilt(:,jj)),exp(otherclusterfilt(:,jj)));
%         pval(jj) = p;
%     end
%     [lowp, sortidx] = sort(pval, 'ascend');
%     sortedgenes = tempfilter(sortidx);
%     disp(uniqueclus(ii));
%     disp(sortedgenes(1:30));
%     sortedfcbyp = log2fc(sortidx);
%     sortedgenes(sortedfcbyp > 4)

writecell(markers,'markers.csv','Delimiter',',');
%% Heatmap
% [sortedcluster sortidx] = sort(cluster, 'ascend');
% sortednormdata = normdata(sortidx,:);
% sortedmarkers = [markers(1:10,:);markers(101:180,:);markers(11:100,:)];
% heatmapdat = sortednormdata(:,ismember(string(filteredgenes),string(sortedmarkers(:,1)')));
% h = HeatMap(flip(heatmapdat'),'RowLabels',sortedmarkers(:,1));
% 
% %% Heatmap
% sortdat = [double(cluster) normdata];
% sortdatbycluster = sortrows(sortdat);
% for ii=1:18
%     heatmapdat(:,(ii-1)*10+1:(ii-1)*10+10) = sortednormdata(:,ismember(string(filteredgenes),string(markers((ii-1)*10+1:(ii-1)*10+10,1)')));
% end

%% Ceramide synthesis gene analysis
% ceramide synthesis pathway genes were determined using Gene Ontology
cergenes = readcell('ceramidesynthesisgenes.txt');
cergenes = [string(cergenes(:,2)); 'FAM57B'];
cerdata = filtereddata(:,ismember(string(filteredgenes),cergenes));

%% Analyze differential gene expression between clusters
n = sum(ismember(string(filteredgenes),cergenes));
cluster = string(idx);
uniqueclus = (string(unique(sort(idx))));
cermarkers = cell(length(uniqueclus)*n, 3);
for ii = 1:length(uniqueclus)

    otherclusters = cerdata(~strcmp(cluster,uniqueclus(ii)),:);
    thiscluster = cerdata(strcmp(cluster,uniqueclus(ii)),:);
    foldchange = mean(exp(thiscluster)) ./ mean(exp(otherclusters));
    log2fc = log2(foldchange);
    [highfc, sortidx] = sort(log2fc, 'descend');
    sortedgenes = cergenes(sortidx);
    cermarkers((ii-1)*n+1:(ii-1)*n+n,1) = cellstr(sortedgenes(1:n)');
    cermarkers((ii-1)*n+1:(ii-1)*n+n,2) = {uniqueclus(ii)};
    cermarkers((ii-1)*n+1:(ii-1)*n+n,3) = num2cell(highfc(1:n));

end
certab = cell2table(cermarkers,'VariableNames',{'Gene','Cluster','logFC'});
h = heatmap(certab,'Cluster','Gene','ColorVariable','logFC');

%% Analyze foveal and peripheral cone ceramide gene differential expression
% extract cone data
conecerdata = cerdata(ismember(string(cluster),'8'),:);
coneregion = region(ismember(string(cluster),'8'),:);
n = sum(ismember(string(filteredgenes),cergenes));
regs = {'fovea' 'peripheral'};
cermarkers = cell(2*n, 4); %initialize result matrix

for ii = 1:2
    otherreg = conecerdata(~strcmp(coneregion,regs(ii)),:);
    thisreg = conecerdata(strcmp(coneregion,regs(ii)),:);
    
    % calculate log fold change
    foldchange = mean(exp(thisreg)) ./ mean(exp(otherreg));
    log2fc = log2(foldchange);
    [highfc, sortidx] = sort(log2fc, 'descend'); % sort by descending fold change
    sortedgenesbyfc = cergenes(sortidx);
    cermarkers((ii-1)*n+1:(ii-1)*n+n,1) = cellstr(sortedgenesbyfc(1:n)');
    cermarkers((ii-1)*n+1:(ii-1)*n+n,2) = {regs(ii)};
    cermarkers((ii-1)*n+1:(ii-1)*n+n,3) = num2cell(highfc(1:n));

    % calculate percent of cells
    otherpercent = sum(otherreg > 0) / size(otherreg,1);
    thispercent = sum(thisreg > 0) / size(thisreg,1);
    deltapercent = thispercent - otherpercent;
    cermarkers((ii-1)*n+1:(ii-1)*n+n,5) = num2cell(deltapercent(sortidx));
    
    % calculate t-test
    pval = zeros(1,size(thisreg,2));
    for jj = 1:size(thisreg,2)
        [h,p] = ttest2(exp(thisreg(:,jj)),exp(otherreg(:,jj)));
        pval(jj) = p;
    end
    cermarkers((ii-1)*n+1:(ii-1)*n+n,4) = num2cell(pval(sortidx));
    disp(cergenes(pval < 0.05)); % displays significant genes
    
end

% Create heatmap
nonzeromarkers = cermarkers(([cermarkers{:,3}]' ~= 0),:); % filter genes with 0 expression
figure;
certab = cell2table(nonzeromarkers,'VariableNames',{'Gene','Region','logFC','pval','percent'});
h = heatmap(certab,'Region','Gene','ColorVariable','logFC');

% Plot logFC by percent of cells
onereg = nonzeromarkers((string([nonzeromarkers{:,2}]') ~= 'peripheral'),:);
figure;
plot([onereg{:,5}]',[onereg{:,3}]','.','MarkerSize',12);
xlabel('delta percent');
ylabel('average logFC');

sigpval = onereg([onereg{:,4}]' < 0.05,:);
text([sigpval{:,5}]'+0.005,[sigpval{:,3}]'+0.03,string(sigpval(:,1)),'FontSize',6);

%% Analyze foveal and peripheral cone gene differential expression for all genes
% extract cone data
conecerdata = filtereddata(ismember(string(cluster),'8'),:);
coneregion = region(ismember(string(cluster),'8'),:);
n = length(string(filteredgenes));
regs = {'fovea' 'peripheral'};
cermarkers = cell(2*n, 4); %initialize result matrix

for ii = 1:2
    otherreg = conecerdata(~strcmp(coneregion,regs(ii)),:);
    thisreg = conecerdata(strcmp(coneregion,regs(ii)),:);
    
    % calculate log fold change
    foldchange = mean(exp(thisreg)) ./ mean(exp(otherreg));
    log2fc = log2(foldchange);
    [highfc, sortidx] = sort(log2fc, 'descend'); % sort by descending fold change
    sortedgenesbyfc = filteredgenes(sortidx);
    cermarkers((ii-1)*n+1:(ii-1)*n+n,1) = cellstr(sortedgenesbyfc(1:n)');
    cermarkers((ii-1)*n+1:(ii-1)*n+n,2) = {regs(ii)};
    cermarkers((ii-1)*n+1:(ii-1)*n+n,3) = num2cell(highfc(1:n));

    % calculate percent of cells
    otherpercent = sum(otherreg > 0) / size(otherreg,1);
    thispercent = sum(thisreg > 0) / size(thisreg,1);
    deltapercent = thispercent - otherpercent;
    cermarkers((ii-1)*n+1:(ii-1)*n+n,5) = num2cell(deltapercent(sortidx));
    
    % calculate t-test
    pval = zeros(1,size(thisreg,2));
    for jj = 1:size(thisreg,2)
        [h,p] = ttest2(exp(thisreg(:,jj)),exp(otherreg(:,jj)));
        pval(jj) = p;
    end
    cermarkers((ii-1)*n+1:(ii-1)*n+n,4) = num2cell(pval(sortidx));
    disp(filteredgenes(pval < 0.05)); % displays significant genes
    
end

% Create heatmap
nonzeromarkers = cermarkers(([cermarkers{:,3}]' ~= 0),:); % filter genes with 0 expression
figure;
certab = cell2table(nonzeromarkers,'VariableNames',{'Gene','Region','logFC','pval','percent'});
h = heatmap(certab,'Region','Gene','ColorVariable','logFC');
%%
% Plot logFC by percent of cells
onereg = nonzeromarkers((string([nonzeromarkers{:,2}]') ~= 'peripheral'),:);
figure;
plot([onereg{:,5}]',[onereg{:,3}]','.','MarkerSize',12);
xlabel('delta percent');
ylabel('average logFC');

sigpval = onereg([onereg{:,3}]' > 3 | [onereg{:,3}]' < -2.5,:);
text([sigpval{:,5}]'+0.005,[sigpval{:,3}]'+0.03,string(sigpval(:,1)),'FontSize',6);