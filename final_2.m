%% Preliminary Analysis
load('meth.mat');
load('microarray.mat');

% Load data individually for ease of access
data_meth = meth.data;
data_gene = microarray.data;
cells_meth = meth.cells;
cells_gene = microarray.cells;
gene_ids_meth = meth.geneIds;
gene_ids_gene = microarray.geneIds;

% Discard first two columns -- not formatted correctly
data_meth = data_meth(:,3:end);
cells_meth = cells_meth(:,3:end);

%% Find genes that are in common between the two sets
[in_common, indexes_meth, indexes_gene] = intersect(gene_ids_meth, gene_ids_gene);
genes_meth = gene_ids_meth(indexes_meth);
data_meth = data_meth(indexes_meth,:);
csvwrite('data_meth.csv', data_meth);
genes_gene = gene_ids_gene(indexes_gene);
data_gene = data_gene(indexes_gene,:);
csvwrite('data_gene.csv', data_gene);

%% Convert gene expression to average gene expression across cell lines
avg_gene = mean(data_gene,2);
%% Linear regression with linear dimensions
avg_meth = mean(data_meth,2);
ones = ones(length(avg_meth),1);
avg_meth_w_ones = [avg_meth ones];
coeffs_basic = avg_meth_w_ones \ avg_gene;
func = @(x) coeffs_basic(1)*x + coeffs_basic(2);
figure(1)
hold on
scatter(avg_meth(:,1), avg_gene);
fplot(func, [0 1], 'LineWidth', 3);
xlabel('Average DNA Methylation Level');
ylabel('Average Gene Expression');
hold off