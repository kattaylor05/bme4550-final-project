%% Preliminary Analysis
load('meth.mat');

% Load data individually for ease of access
data = meth.data;
cells = meth.cells;
gene_ids = meth.geneIds;

% Discard first two columns -- not formatted correctly
data = data(:,3:end);
cells = cells(:,3:end);

% Find average methylation for each gene
avg_meth = mean(data,2);
figure(1)
hist(avg_meth);
xlabel("Average Methylation");
ylabel("Occurrences");
[sorted_avg, indexes] = sort(avg_meth);
genes_sorted_mean = gene_ids(indexes);
data_sorted_mean = data(indexes,:);
genes_high_meth = genes_sorted_mean(24951:end);
data_high_meth = data_sorted_mean(24951:end,:);
genes_low_meth = genes_sorted_mean(1:4244);
data_low_meth = data_sorted_mean(1:4244,:);

% Find methylation variability for each gene
var_meth = var(data,0,2);
figure(2)
hist(var_meth);
xlabel("Variance in Methylation");
ylabel("Occurrences");
[sorted_var, indexes_var] = sort(var_meth);
genes_sorted_var = gene_ids(indexes_var);
data_sorted_var = data(indexes,:);
genes_high_var = genes_sorted_var(25666:end);
data_high_var = data_sorted_var(25666:end,:);

%% Find pathways
load('p53Data_HW2.mat');
[in_common, index_ids, index_names] = intersect(gene_ids, s.gene.names);
path_data = s.msig.data(index_names,:);
common_data = data(index_ids,:);

%% Low Methylation
[in_common_low, index_a_low, index_b_low] = intersect(genes_low_meth, in_common, 'stable');
path_data_low = path_data(index_b_low,:);
common_data_low = data_low_meth(index_a_low,:);
p_vals = zeros(522,1);
for i = 1:522
    pathway_hits = path_data_low(:,i);
    index_hits = find(pathway_hits);
    [es_scores, max_es, index_low] = es_score(mean(common_data_low,2), index_hits, 1:1847, 1847, 1);
    max_nulls = zeros(100,1);
    for j = 1:100
        random_hits = randi(1847, length(index_hits), 1);
        random_hits = sort(random_hits);
        [es_scores_null, max_es_null, index_null] = es_score(mean(common_data_low,2), random_hits, 1:1847, 1847, 1);
        max_nulls(j) = max_es_null;
    end
    test_stat = (max_es - mean(max_nulls)) / std(max_nulls);
    p = 1 - normcdf(test_stat);
    p_vals(i) = p;
end
figure(3)
hist(max_nulls)
%% Processing Low Methylated Significant Pathways
[sorted_pvals, indexes_pvals] = sort(p_vals);
pathways_low = (s.msig.names(indexes_pvals(1:25)))';
disp(pathways_low)
%% High Methylation
[in_common_high, index_a_high, index_b_high] = intersect(genes_high_meth, in_common);
path_data_high = path_data(index_b_high,:);
common_data_high = data_high_meth(index_a_high,:);
p_vals_high = zeros(522,1);
for i = 1:522
    pathway_hits = path_data_high(:,i);
    index_hits = find(pathway_hits);
    [es_scores, max_es, index_low] = es_score(mean(common_data_high,2), index_hits, 1:1004, 1004, 1);
    max_nulls = zeros(100,1);
    for j = 1:100
        random_hits = randi(1004, length(index_hits), 1);
        random_hits = sort(random_hits);
        [es_scores_null, max_es_null, index_null] = es_score(mean(common_data_high,2), random_hits, 1:1004, 1004, 1);
        max_nulls(j) = max_es_null;
    end
    test_stat = (max_es - mean(max_nulls)) / std(max_nulls);
    p = 1 - normcdf(test_stat);
    p_vals_high(i) = p;
end
%% Processing Data for High Methylation
[sorted_pvals_high, indexes_pvals_high] = sort(p_vals_high);
pathways_high = (s.msig.names(indexes_pvals_high(1:54)))';
disp(pathways_high);
%% High Variance
[in_common_var, index_a_var, index_b_var] = intersect(genes_high_var, in_common);
path_data_var = path_data(index_b_var,:);
common_data_var = data_high_var(index_a_var,:);
p_vals_var = zeros(522,1);
for i = 1:522
    pathway_hits = path_data_var(:,i);
    index_hits = find(pathway_hits);
    [es_scores, max_es, index_low] = es_score(mean(common_data_var,2), index_hits, 1:762, 762, 1);
    max_nulls = zeros(100,1);
    for j = 1:100
        random_hits = randi(1004, length(index_hits), 1);
        random_hits = sort(random_hits);
        [es_scores_null, max_es_null, index_null] = es_score(mean(common_data_var,2), random_hits, 1:762, 762, 1);
        max_nulls(j) = max_es_null;
    end
    test_stat = (max_es - mean(max_nulls)) / std(max_nulls);
    p = 1 - normcdf(test_stat);
    p_vals_var(i) = p;
end
%% Processing Data for High Methylation
[sorted_pvals_var, indexes_pvals_var] = sort(p_vals_var);
pathways_var = (s.msig.names(indexes_pvals_var(1:8)))';