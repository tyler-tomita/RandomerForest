function T = randtile(ImHeight,ImWidth,MaxHeight,MaxWidth,N)
% RANDTILE   Sample random tile of contiguous pixels from an image
%
%   T = RANDTILE returns an ImHeight x ImWidth x N array of logical zeros 
%   and ones, where ones comprise a contiguous randomly sized and placed
%   tile in each of the N slices. MaxHeight and MaxWidth control the
%   maximum height and width of the tiles respectively

if isempty(MaxHeight)
    MaxHeight = ImHeight;
end

if isempty(MaxWidth)
    MaxWidth = ImWidth;
end

T = false(ImHeight,ImWidth);

for i = 1:N
    % sample random height and width of the tile
    TileHeight = randi(MaxHeight);
    TileWidth = randi(MaxWidth);
    
    % sample random coordinates for the top-left of the tile
    RowIdx = randi(ImHeight-TileHeight+1);
    ColIdx = randi(ImWidth-TileWidth+1);
    
    T(RowIdx:RowIdx+TileHeight-1,ColIdx:ColIdx+TileWidth-1,i) = true;
end