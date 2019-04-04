function data_zscore = zscore_4D(data);
% Compute zscores of 4D matrix

dim = size(data);
if length(dim) ~= 4; error('Input must to be a 4D matrix'); end
data_zscore = zeros(dim(1),dim(2),dim(3),dim(4));

for i = 1:size(data, 1)
    for j = 1:size(data, 2)
        for k = 1:size(data, 3)
            a = reshape(data(i,j,k,:), 1, size(data, 4));
            if std(a)>0
                data_zscore(i,j,k,:) = zscore(a);
            end
        end
    end
end

end