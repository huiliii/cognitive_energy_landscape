function build_process_profiles(expression_file, gene_list_file, output_file)
% Construct six process-associated AHBA profiles using the study workflow.
% Inputs:
%   expression_file: preprocessed AAL2 region-by-gene CSV (label + genes).
%   gene_list_file: process,gene CSV supplied alongside this script.
%   output_file: destination MAT file.
% Use the SAME region and gene ordering as the structural connectomes.

AHBA = readtable(expression_file, 'VariableNamingRule', 'preserve');
gene_lists = readtable(gene_list_file, 'TextType', 'string');
process_names = ["Synapse development", "Neuron migration", ...
    "Neuron differentiation", "Dendrite development", ...
    "Axon development", "Myelination"];

% Original AAL2 row selection (MATLAB indexing). The original matrix has
% 120 regions, with its label column ordered 1..120. Do not apply these
% positions to another atlas or a reordered matrix without checking labels.
assert(height(AHBA) == 120, 'Expected the original 120-row AAL2 matrix.');
assert(strcmp(AHBA.Properties.VariableNames{1}, 'label'), ...
    'The first column must contain AAL2 labels.');
assert(isequal(AHBA.label(:), (1:120)'), 'Expected AAL2 labels 1..120.');
aal_index = [1:28 33:94];
region_labels = AHBA.label(aal_index);
AHBA_gene_names = AHBA.Properties.VariableNames(2:end)';
all_neuro_events = cell(7, 1);
for p = 1:6
    all_neuro_events{p+1} = cellstr(unique(gene_lists.gene( ...
        gene_lists.process == process_names(p))));
end
all_neuro_events{1} = unique(vertcat(all_neuro_events{2:end}));

brain_nd_rna_map = cell(7, 1);
matched_genes = cell(7, 1);
unmatched_genes = cell(7, 1);
for p = 1:7
    candidates = unique(all_neuro_events{p});
    [matched_genes{p}, gene_idx] = intersect(AHBA_gene_names, candidates);
    unmatched_genes{p} = setdiff(candidates, AHBA_gene_names);
    assert(~isempty(gene_idx), 'A process has no matched genes.');
    brain_nd_rna_map{p} = table2array(AHBA(aal_index, gene_idx+1));
    fprintf('Profile %d: %d matched genes, %d unmatched genes\n', ...
        p, numel(matched_genes{p}), numel(unmatched_genes{p}));
end

% Missing-region handling must be explicit. The historical intermediate
% matrix has missing expression at AAL2 label 92, while the saved final
% profile there matches label 91. That intermediate edit is not recorded
% in the extraction script, so it is not silently applied here.
assert(all(cellfun(@(x) all(isfinite(x), 'all'), brain_nd_rna_map)), ...
    ['Missing expression detected. Document and apply your chosen missing-' ...
     'region handling before constructing process profiles.']);
brain_nd_rna_map_mean = cellfun(@(x) mean(x, 2), ...
    brain_nd_rna_map, 'UniformOutput', false);
brain_nd_rna_map_mean = cat(2, brain_nd_rna_map_mean{:});
% Column 1 averages the unique union of matched genes. All profiles share
% the same scale, preserving between-process differences in expression.
common_scale = mean(brain_nd_rna_map_mean(:, 1), 'all');
brain_nd_rna_map_nm = brain_nd_rna_map_mean ./ common_scale;
profile_names = ["Union of six gene sets", process_names];
save(output_file, 'brain_nd_rna_map', 'brain_nd_rna_map_mean', ...
    'brain_nd_rna_map_nm', 'region_labels', 'profile_names', ...
    'matched_genes', 'unmatched_genes', 'common_scale');
end
