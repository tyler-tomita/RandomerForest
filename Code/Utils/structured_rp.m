function [A,ClassPairIdx] = structured_rp(ImHeight,ImWidth,MaxHeight,MaxWidth,d,Mu_diff,Control)
% STRUCTURED_RP   Random projections with nonzero entries corresponding to
% contiguous rectangular tiles of pixels
%
%   A = STRUCTURED_RP(IMHEIGHT,IMWIDTH,MAXHEIGHT,MAXWIDTH,D,MU_DIFF,CONTROL)
%   returns an ImHeight*ImWidth x d random matrix with nonzero entries in
%   each column restricted to a contiguous tile of pixels. MaxHeight and
%   MaxWidth specify the maximum allowable height and width of image tiles.
%   Mu_diff is a p x K-1 matrix in which columns are the mean difference
%   vectors of each of the K-1 pairs of classes. If Mu_diff is supplied,
%   then they are used in forming the structured random matrix A.
%   Otherwise, A is formed by sampling -1 w/ prob 0.5 and +1 w/ prob 0.5
%   for each nonzero entry. Control is a logical indicating if the control
%   version of random projections should be executed. The
%   control version generates a random matrix of the same density as
%   structured random projections, but does not spatially constrain the
%   construction of A.

p = ImHeight*ImWidth;
A = zeros(p,d);

if isempty(MaxHeight)
    MaxHeight = ImHeight;
end

if isempty(MaxWidth)
    MaxWidth = ImWidth;
end

isdiff = ~isempty(Mu_diff);

npairs = size(Mu_diff,2);

if isdiff
    if npairs > 1
        ClassPairIdx = randi(npairs,1,d);
    else
        ClassPairIdx = ones(1,d);
    end
else
    ClassPairIdx = [];
end

% sample random heights and widths of the tiles
TileHeight = randi(MaxHeight,1,d);
TileWidth = randi(MaxWidth,1,d);
if ~Control
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
            A(NzIdx(:),i) = Mu_diff(NzIdx(:),ClassPairIdx(i));
        end
    end
else
    nnzs = TileHeight.*TileWidth;
    for i = 1:d
        NzIdx = randperm(p,nnzs(i));
        A(NzIdx,i) = round(rand(length(NzIdx),1))*2-1;
    end
end

A = sparse(A);