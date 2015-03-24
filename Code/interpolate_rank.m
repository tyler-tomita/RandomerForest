function Xnew_rank = interpolate_rank(X,Xnew)
    Xnew_rank = zeros(size(Xnew));
    Xrank = passtorank(X);
    for i = 1:size(Xnew,1)
        Xcat = cat(1,X,Xnew(i,:));
        Xcat_rank = passtorank(Xcat);
        for j = 1:size(Xnew(i,:),2)
            if Xcat_rank(end,j) == size(Xcat_rank,1) || Xcat_rank(end,j) == 1
                Xnew_rank(i,j) = Xcat_rank(end,j);
            else
                lower_idx = find(Xrank(:,j)==Xcat_rank(end,j)-1);
                upper_idx = find(Xrank(:,j)==Xcat_rank(end,j)+1);
                Xnew_rank(i,j) = ((Xnew(i,j) - X(lower_idx,j))/(X(upper_idx,j) - X(lower_idx,j)))*(Xrank(upper_idx,j) - Xrank(lower_idx,j)) + Xrank(lower_idx,j);
            end
        end
    end
end