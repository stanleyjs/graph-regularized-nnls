function [Y,library_sizes] = library_size_normalization(X,median_scaled)
% implements a library size normalization of X
% l1 normalize columns of input, then multiply by the median l1 norm if
% `median_scaled`

if nargin<2
    median_scaled = true;
end
library_sizes = sum(abs(X),1);
Y = X./library_sizes;
if median_scaled
    Y = Y.*median(library_sizes);
end
end

