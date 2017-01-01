function [R,S,SortIdx] = passtorank(A)
%Passes each column of a matrix A to ranks

%A: n x d matrix
%R: represents the matrix obtained by converting each element of A to its
%column marginal rank
    [S,SortIdx] = sort(A,1);
    [nrow,ncol] = size(SortIdx);
    i_adjust = repmat(0:nrow:(ncol-1)*nrow,nrow,1);
    SortIdx_adj = SortIdx + i_adjust;
    R = A;
    R(SortIdx_adj) = repmat(1:size(R,1),1,size(R,2));
end