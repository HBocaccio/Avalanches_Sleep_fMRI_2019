function cluster_sizes = get_cluster_sizes(clusters);
% Get clusters sizes from 4D matrix of clusters labels
%
% Inputs
%    clusters - matrix of clusters labels, allowing formats:
%               4D matrix ('clusters_labeling' output)
%               Nx5 sparse matrix ('clusters2sparse' output)
%               1 x time cell matrix ('sparse2cell' output)
%
% Outputs
%    clusters_sizes - Variable containing clusters sizes with columns [size, label, time]
%
% Usage: 
%    after 'clusters_labeling'
%
% Reference:
% Tagliazucchi, E., Balenzuela, P., Fraiman, D., & Chialvo, D. R. (2012). 
% Criticality in large-scale brain fMRI dynamics unveiled by a novel 
% point process analysis. Frontiers in physiology, 3, 15.
%

dim = size(clusters);
if length(dim) == 2 && ~iscell(clusters); 
    representation = 'sparse';
elseif length(dim) == 2 && iscell(clusters); 
    representation = 'cell';
elseif length(dim) > 2; 
    representation = '4d'; 
else error('Input Dimension Error'); 
end

switch representation
    case '4d'
        data_bis = reshape(clusters, dim(1)*dim(2)*dim(3), dim(4));
        [ids_cell, ~, ids_positions_cell] = arrayfun(@(k) unique(data_bis(data_bis(:, k) ~= 0, k), 'stable'), 1:size(data_bis, 2), 'un', 0);
        sizes_cell = arrayfun(@(k) accumarray(ids_positions_cell{k}, data_bis(data_bis(:, k) ~= 0, k), size(ids_cell{k}), @(x) numel(x)), 1:size(data_bis, 2), 'un', 0);
    case 'sparse'
        data_bis = [clusters(:, 1), clusters(:, 5)];
        [ids_cell, ~, ids_positions_cell] = arrayfun(@(k) unique(data_bis(data_bis(:, 2) == k, 1), 'stable'), 1:max(data_bis(:,2)), 'un', 0);        
        sizes_cell = arrayfun(@(k) accumarray(ids_positions_cell{k}, data_bis(data_bis(:, 2) == k, 1), size(ids_cell{k}), @(x) numel(x)), 1:max(data_bis(:,2)), 'un', 0);
    case 'cell'
        [ids_cell, ~, ids_positions_cell] = arrayfun(@(k) unique(clusters{k}(:,1), 'stable'), 1:size(clusters,2), 'un', 0);
        sizes_cell = arrayfun(@(k) accumarray(ids_positions_cell{k}, clusters{k}(:,1), size(ids_cell{k}), @(x) numel(x)), 1:size(clusters,2), 'un', 0);
end

csizes_temp_cell = cellfun(@sortrows, cellfun(@horzcat, sizes_cell, ids_cell, 'UniformOutput',false), 'UniformOutput',false);
n_ids = cellfun(@numel, ids_cell);
t_cell = arrayfun(@(k) repmat(k, n_ids(k), 1), 1:length(n_ids), 'un', 0);
cluster_sizes = [vertcat(csizes_temp_cell{:}), vertcat(t_cell{:})];

end