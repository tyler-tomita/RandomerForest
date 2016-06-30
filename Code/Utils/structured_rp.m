function A = structured_rp(ImHeight,ImWidth,MaxHeight,MaxWidth,d,mu_diff)
% STRUCTURED_RP   Random projections with nonzero entries corresponding to
% contiguous rectangular tiles of pixels
%
%   A = STRUCTURED_RP returns an ImHeight*ImWidth x d random matrix with
%   nonzero entries in each column restricted to a contiguous tile of
%   pixels

p = ImHeight*ImWidth;
A = zeros(p,d);

if isempty(MaxHeight)
    MaxHeight = ImHeight;
end

if isempty(MaxWidth)
    MaxWidth = ImWidth;
end

isdiff = ~isempty(mu_diff);

npairs = size(mu_diff,2);

if isdiff
    if npairs > 1
        ClassPairIdx = randi(npairs,1,d);
    else
        ClassPairIdx = ones(1,d);
    end
end

% sample random heights and widths of the tiles
TileHeight = randi(MaxHeight,1,d);
TileWidth = randi(MaxWidth,1,d);

if ~isdiff
    for i = 1:d
        % sample random coordinates for the top-left of the tile
        PossibleIdx = repmat((1:ImHeight-TileHeight(i)+1)',1,...
            ImWidth-TileWidth(i)+1) +...
            repmat(0:ImHeight:ImHeight*(ImWidth-TileWidth(i)),...
            ImHeight-TileHeight(i)+1,1);
        Element = randi(numel(PossibleIdx));
        TopLeftIdx = PossibleIdx(Element);

        %Sample -1 or +1 for each nonzero coordinate
        NzIdx = repmat((TopLeftIdx:TopLeftIdx+TileHeight(i)-1)',1,TileWidth(i))...
            + repmat(0:ImHeight:ImHeight*(TileWidth(i)-1),TileHeight(i),1);
        A(NzIdx(:),i) = round(rand(numel(NzIdx),1))*2-1;
    end
else
    for i = 1:d
        % sample random coordinates for the top-left of the tile
        PossibleIdx = repmat((1:ImHeight-TileHeight(i)+1)',1,...
            ImWidth-TileWidth(i)+1) +...
            repmat(0:ImHeight:ImHeight*(ImWidth-TileWidth(i)),...
            ImHeight-TileHeight(i)+1,1);
        Element = randi(numel(PossibleIdx));
        TopLeftIdx = PossibleIdx(Element);

        %Replace nonzero coordinates with mean difference of a random pair
        %of class labels
        NzIdx = repmat((TopLeftIdx:TopLeftIdx+TileHeight(i)-1)',1,TileWidth(i))...
            + repmat(0:ImHeight:ImHeight*(TileWidth(i)-1),TileHeight(i),1);
        A(NzIdx(:),i) = mu_diff(NzIdx(:),ClassPairIdx(i));
    end
end

A = sparse(A);