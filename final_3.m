%% Preliminary Analysis
load('meth.mat');
load('microarray.mat');
load('RPPA.mat');

% Load data individually for ease of access
data_meth = meth.data;
data_gene = microarray.data;
data_prot = RPPA.data;
cells_meth = meth.cells;
cells_gene = microarray.cells;
cells_prot = RPPA.cells;
gene_ids_meth = meth.geneIds;
gene_ids_gene = microarray.geneIds;
proteins = RPPA.proteins;

% Discard first two columns -- not formatted correctly
data_meth = data_meth(:,3:end);
cells_meth = cells_meth(:,3:end);

% Discard cell line with no invasiveness data for proteins
to_delete = [1 2 11 12 15];
data_prot(:,to_delete) = [];
cells_prot(:,to_delete) = [];

% Discard cell lines with no invasiveness for gene expression
to_delete_2 = [1 13 14 16];
data_gene(:,to_delete_2) = [];
cells_gene(:,to_delete_2) = [];

% Invasiveness data
invasive_meth = [0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 0 1 0 1 0 0 0 0 0 1 1 0 0 1 0 0 0 0];
invasive_gene = [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 1 1 0 0 1 0 1 0 0 0 0 0 1 1 0 0 1 0 0 0 0 0 1 1];
invasive_prot = [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 1 0 1 1 1 0 1 1 0 0 0 0 1 1 1 0 0 0 0 0];

% Write to csv for Python analysis
csvwrite('data_meth_full.csv', data_meth);
csvwrite('data_gene_full.csv', data_gene);
csvwrite('data_protein.csv', data_prot);
csvwrite('invasive_meth.csv', invasive_meth);
csvwrite('invasive_gene.csv', invasive_gene);
csvwrite('invasive_prot.csv', invasive_prot);