function pvalue = permutation_test(x, y, nperm)
% Permutation test to compare means of distributions

if nargin < 3 || isempty(nperm)
    nperm = 1000;
end

x = x(~isnan(x));
y = y(~isnan(y));
x = x(:);
y = y(:);
z = [x; y];
n1 = length(x);
n2 = length(y);

statistic = mean(x) - mean(y);
[~, idx] = sort(rand(nperm, n1+n2), 2);
x_perm = z(idx(:, 1:n1));
y_perm = z(idx(:, n1+1:end));
statistic_perm = mean(x_perm, 2) - mean(y_perm, 2);

b = sum(abs(statistic_perm) >= abs(statistic));
pvalue = (b + 1) / (nperm + 1);

end