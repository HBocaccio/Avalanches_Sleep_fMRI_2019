function std_data = std_4D(data)
% Compute std of 4D matrix
for i = 1:size(data, 1)
    for j = 1:size(data, 2)
        for k = 1:size(data, 3)
            a = reshape(data(i, j, k, :), 1, size(data, 4));
            std_data(i, j, k) = std(a);
        end
    end
end
end