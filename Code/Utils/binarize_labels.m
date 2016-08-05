function Yb = binarize_labels(Y,Labels)
% BINARIZE(Y) encodes an N x 1 array of multi-class labels Y into an N x
% n_classes binary matrix Yb. Labels is the set of unique class labels. The
% columns of Yb reflect the order of Labels.

n_classes = length(Labels);

Yb = zeros(length(Y),n_classes);

if iscell(Y) && iscell(Labels)
    for i = 1:n_classes
        Yb(strcmp(Y,Labels{i}),i) = 1;
    end
else
    for i = 1:n_classes
        Yb(Y==Labels(i),i) = 1;
    end
end

end