function sparse_clusters_cell = sparse2cell(sparse_clusters)
% Transform Nx5 2D matrix sparse representation to cell (1 x time)

dim = size(sparse_clusters);
if length(dim) ~= 2; error('Input must to be a 2D sparse matrix'); end

sparse_clusters_cell = arrayfun(@(k) sparse_clusters(sparse_clusters(:,5) == k, 1:4), 1:max(sparse_clusters(:,5)), 'un', 0);

end