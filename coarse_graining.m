function data_downsampled = coarse_graining(data);
% Obtain spatial downsampling of the data through a coarse graining method
%
% Inputs
%    data - 4D matrix with time series (x,y,z,t)
%
% Outputs
%    data_downsampled - 4D matrix with downsampled data
%
% Reference:

dim = size(data);
if length(dim) ~= 4; error('Input must to be a 4D matrix'); end

for t = 1:size(data, 4)
    fprintf('Progress: %f %%\n' ,t/size(data, 4)*100)
    volume = data(:, :, :, t);
    volume_subsampled = zeros(floor(size(volume, 1)/2), floor(size(volume, 2)/2), floor(size(volume, 3)/2));
    for i = 1:floor(size(volume,1)/2)
        for j = 1:floor(size(volume,2)/2)
            for k = 1:floor(size(volume,3)/2)
                block = volume(2*i-1:2*i, 2*j-1:2*j,2*k-1:2*k);
                volume_subsampled(i, j, k) = mean(block(:));
            end
        end
    end
    data_downsampled(:, :, :, t) = volume_subsampled;
end

end