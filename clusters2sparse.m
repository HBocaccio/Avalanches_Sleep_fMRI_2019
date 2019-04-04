function sparse_clusters = clusters2sparse(clusters);
% Transform 4D matrix of clusters labels to sparse representation

dim = size(clusters);
if length(dim) ~= 4; error('Clusters matrix dimension error'); end

sparse_clusters = zeros(length(find(clusters)), 5);
e = 1;
for t = 1:dim(4)
    for zid = 1:dim(3)
        [x, y, v] = find(clusters(:, :, zid, t));
        if length(v) ~= 0
            sparse_clusters(e:e + length(v) - 1, :) = [v, x, y, repmat(zid, length(v), 1), repmat(t, length(v), 1)];
            e = e + length(v);
        end
    end    
end
    
end