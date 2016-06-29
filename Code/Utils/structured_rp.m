function A = structured_rp(ImHeight,ImWidth,MaxHeight,MaxWidth,d)
% STRUCTURED_RP   Random projections with nonzero entries corresponding to
% contiguous rectangular tiles of pixels
%
%   A = STRUCTURED_RP returns an ImHeight*ImWidth x d random matrix with
%   nonzero entries in each column restricted to a contiguous tile of
%   pixels

p = ImHeight*ImWidth;
A = sparse(p,d);

if isempty(MaxHeight)
    MaxHeight = ImHeight;
end

if isempty(MaxWidth)
    MaxWidth = ImWidth;
end

for i = 1:d
    % sample random height and width of the tile
    TileHeight = randi(MaxHeight);
    TileWidth = randi(MaxWidth);
    
    % sample random coordinates for the top-left of the tile
    RowIdx = randi(ImHeight-TileHeight+1);
    ColIdx = randi(ImWidth-TileWidth+1);
    [r,c] = meshgrid(RowIdx:RowIdx+TileHeight-1,ColIdx:ColIdx+TileWidth-1);
    nnzs = sub2ind([ImHeight,ImWidth],r(:),c(:));
    A(nnzs,i) = round(rand(length(nnzs),1))*2-1;
end