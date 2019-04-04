function data_shuffled = phase_shuffling_4D(data, mask);
% Add a random phase to time series from 4D matrix
%
% Inputs
%    data - 4D matrix with time series (x,y,z,t)
%    mask (optional) - Mask with same spatial dimensions and normalization
%                      (e.g. MNI) than 'data'. Default mask obtained from
%                      voxels of time series with std == 0 (constant).
%
% Outputs
%    data_shuffled - 4D matrix with surrogate time series
%
% Reference:

dim = size(data);
if length(dim) ~= 4; error('Input must to be a 4D matrix'); end

if nargin == 1
    mask = std_4D(data);
end

data_shuffled = zeros(size(data,1),size(data,2),size(data,3),size(data,4)-1);
e = 1;
for x = 1:size(data, 1)
    for y = 1:size(data, 2)
        for z = 1:size(data, 3)
            if mask(x, y, z) ~= 0
                fprintf('Progress: %f %%\n' ,e/numel(find(mask ~= 0))*100)
                ts = reshape(data(x, y, z, :), 1, size(data, 4));
                ts_shuffled = phase_shuffling(ts);
                data_shuffled(x, y, z, :) = ts_shuffled;
                e = e + 1;
            end
        end
    end
end

end

function std_data = std_4D(data)
for i = 1:size(data, 1)
    for j = 1:size(data, 2)
        for k = 1:size(data, 3)
            a = reshape(data(i, j, k, :), 1, size(data, 4));
            std_data(i, j, k) = std(a);
        end
    end
end
end