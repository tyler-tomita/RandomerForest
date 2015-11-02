function R = passtorank(A)
%Passes each column of a matrix A to ranks

%A: n x d matrix
%R: represents the matrix obtained by converting each element of A to its
%column marginal rank
    [S,i] = sort(A,1);
    i_adjust = repmat(0:size(i,1):(size(i,2)-1)*size(i,1),size(i,1),1);
    i = i + i_adjust;
    R = A;
    R(i) = repmat(1:size(R,1),1,size(R,2));
end