# Avalanches_Sleep_fMRI_2019

### Main functions:
- 'clusters_labeling', detect and label clusters of contiguous co-activated voxels (i.e. super-threshold voxels according to point-process transformation).
- 'clusters2sparse', transform 4D matrix of clusters labels to sparse representation.
- 'sparse2cell', transform Nx5 2D matrix sparse representation to cell (1 x time).
- 'get_cluster_sizes', get clusters sizes from 4D matrix of clusters labels.
- 'coarse_graining', obtain spatial downsampling of the data through a coarse graining method.
- 'phase_shuffling_4D', add a random phase to time series from 4D matrix.
- 'permutation_test', permutation test to compare means of distributions.

### Secondary functions:
- 'std_4D', compute std of 4D matrix.
- 'zscore_4D', compute zscores of 4D matrix.
- 'get_position', get single-index target from 3D position.
- 'phase_shuffling', add a random phase to single time serie.

### Workflow:

From 4D matrix of fMRI data containing all BOLD series, we first detected and labeled the clusters of co-activated voxels by applying a point-process transformation and computing connected components in a first-neighbors network, using 'clusters_labeling' function. The output variable is a matrix with the same dimensions than data, but with the respective clusters labels in the position of super-threshold voxels for each volume (instead of BOLD values) and zero values for inactive voxels. This matrix could be transformed to a sparse matrix representation through 'clusters2sparse' function (and then to a cell variable from the sparse matrix of each volume using 'sparse2cell'). Both original and sparse representations (either matrix or cell variables) could be the inputs of 'get_cluster_sizes' function. We used this function to obtain the cluster sizes distributions.

We then grouped sizes for different sleep stages and applied Clauset et al algorithm to each distribution, estimating power-law fit parameters.

The statistical comparisons were made mostly using MATLAB stats toolbox functions ('kruskalwallis', 'ranksum', 'anova1'), exceptuating the permutation tests that were made using 'permutation_test' function.

The same analysis was applied to the spatial down-sampled data, obtained by progressively averaging the signals in voxel blocks, i.e. using a coarse graining technique in volume. Down-sampled data was obtained through the 'coarse_graining' function.

BOLD signals were phase shuffled to construct null models using signals that contain the same spectrum as the original versions, but with scrambled phases (surrogates). This procedure was implemented by applying a Fast Fourier Transform (FFT) to transform real BOLD signals into the frequency domain, and subsequently reversing the FFT after adding a random phase to obtain the real surrogate time series. Phase-shuffled data was obtained using the 'phase_shuffling_4D' function from 4D matrix as input variable. 'phase_shuffling_4D' adds random phases to time series one-by-one, calling the 'phase_shuffling' function to obtain each surrogate time serie.

### Comments:

The 'clusters_labeling' function uses nested functions to compute the standard deviation of 4D data ('std_4D'), z-score normalization of 4D data ('zscore_4D'), and single-index target from 3D position ('get_position'). This functions were also included in the repository as secondary functions, allowing to apply them beyond clusters detection and labeling, and also making easy the manually run of function by sections.

### How to cite:

The point-process transformation in fMRI data was first described and developed in
- Tagliazucchi, E., Balenzuela, P., Fraiman, D., & Chialvo, D. R. (2012). Criticality in large-scale brain fMRI dynamics unveiled by a novel point process analysis. Frontiers in physiology, 3, 15.

The codes and functions uploaded in this repository are an update to apply in https://github.com/HBocaccio/Avalanches_Sleep_fMRI_2019
