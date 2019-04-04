function clusters = clusters_labeling(data, threshold, mask);
% Detect and label clusters of contiguous co-activated voxels 
% (i.e. super-threshold voxels according a point-process transformation)
%
% Requires matlab bgl library for the connected components computation
% (https://github.com/dgleich/matlab-bgl)
%
% Inputs
%    data - 4D matrix (3 spatial dimensions + time). Pre-processed 4D fMRI image   
%    threshold (optional) - Point-process threshold value (s.d.). Default
%                           value of threshold set to 1 (s.d.).
%    mask (optional) - Mask with same spatial dimensions and normalization
%                      (e.g. MNI) than 'data'. Default mask obtained from
%                      voxels of time series with std == 0 (constant).
%
% Outputs
%    clusters - 4D matrix of clusters labels 
%
% Reference:
% Tagliazucchi, E., Balenzuela, P., Fraiman, D., & Chialvo, D. R. (2012). 
% Criticality in large-scale brain fMRI dynamics unveiled by a novel 
% point process analysis. Frontiers in physiology, 3, 15.
%

dim = size(data);
if length(dim) ~= 4; error('Input must to be a 4D matrix'); end

if ~exist('threshold','var'); threshold = 1; end
if ~exist('mask','var'); mask = std_4D(data); end

xmax = dim(1); ymax = dim(2); zmax = dim(3); time = dim(4);
clusters = zeros(xmax,ymax,zmax,time);
Bt = [];

disp('getting dictionary...');
dict = zeros(xmax*ymax*zmax,3);
count = 0;
for x = 1:xmax
    for y = 1:ymax
        for z = 1:zmax
            count = count + 1;
            dict(count,1) = x;
            dict(count,2) = y;
            dict(count,3) = z;
        end
    end
end

% Compute zscores
a = zscore_4D(data);

disp('deleting voxels outside of mask...')
for x = 1:xmax
    for y = 1:ymax
        for z = 1:zmax
            if mask(x, y, z) == 0
                a(x, y, z, :) = zeros(1, time);
            end
        end
    end
end

for t = 1:1:time
    t
    tic
    disp('making graph...')
    pos1 = 0; pos2 = 0; pos3 = 0; x1 = []; y1 = [];
    % 'posx' is used for isolated voxels, assigning labels from 30000
    posx = 30000;

    tic
    for i = 2:xmax-1
        for j = 2:ymax-1
            for k = 2:zmax-1
                n = get_position(i,j,k,xmax,ymax,zmax); 
                if a(i,j,k,t) > threshold && a(i+1,j,k,t) < threshold && a(i-1,j,k,t) < threshold && a(i,j+1,k,t) < threshold && a(i,j-1,k,t) < threshold && a(i,j,k+1,t) < threshold && a(i,j,k-1,t) < threshold
                    clusters(i,j,k,t) = posx;
                    posx = posx+1;
                end
                if a(i,j,k,t) > threshold
                    if a(i+1,j,k,t) > threshold
                        pos1 = pos1+1;
                        x1(pos1) = n;
                        y1(pos1) = get_position(i+1,j,k,xmax,ymax,zmax);
                    end
                    if a(i-1,j,k,t) > threshold 
                        pos1 = pos1+1;
                        x1(pos1) = n;
                        y1(pos1) = get_position(i-1,j,k,xmax,ymax,zmax);
                    end
                    if a(i,j+1,k,t) > threshold
                        pos1 = pos1+1;
                        x1(pos1) = n;
                        y1(pos1) = get_position(i,j+1,k,xmax,ymax,zmax);
                    end
                    if a(i,j-1,k,t) > threshold
                        pos1 = pos1+1;
                        x1(pos1) = n;
                        y1(pos1) = get_position(i,j-1,k,xmax,ymax,zmax);
                    end
                    if a(i,j,k+1,t) > threshold
                        pos1 = pos1+1;
                        x1(pos1) = n;
                        y1(pos1) = get_position(i,j,k+1,xmax,ymax,zmax);
                    end
                    if a(i,j,k-1,t) > threshold
                        pos1 = pos1+1;
                     	x1(pos1) = n;
                        y1(pos1) = get_position(i,j,k-1,xmax,ymax,zmax);
                    end
                end
            end
        end
    end
    toc
    
    if not(isempty(x1)) && not(isempty(y1)) 
        S1 = sparse(x1,y1,ones(1,length(x1)), xmax*ymax*zmax,xmax*ymax*zmax);
        disp('detecting connected components...');
        B1 = components(S1);

        % To avoid variable overwriting
        Bl(t) = length(B1);
        Bt = [Bt B1];

        count = 0;
        for i = 1:1:max(B1)
            l = find(B1 == i) ;
            if length(l) > 1
                count = count+1;
                for i = 1:length(l)
                    clusters(dict(l(i),1),dict(l(i),2),dict(l(i),3),t) = count;
                end
            end
        end
    end
    toc
end
end

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

function data_zscore = zscore_4D(data);
% Compute zscores of 4D matrix
dims = size(data);
data_zscore = zeros(dims(1),dims(2),dims(3),dims(4));
for ii = 1:dims(1)
    for jj = 1:dims(2)
        for kk = 1:dims(3)
            aa = reshape(data(ii,jj,kk,:), 1, dims(4));
            if std(aa)>0
                data_zscore(ii,jj,kk,:) = zscore(aa);
            end
        end
    end
end
end

function output = get_position(i,j,k,xmax,ymax,zmax)
% Get single-index target from 3D position
output = ymax*zmax*(i-1) + zmax*(j-1) + k;
end   