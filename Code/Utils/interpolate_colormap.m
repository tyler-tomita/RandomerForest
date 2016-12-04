function C = interpolate_colormap(ColorMap,nhues,TwoColor)
    %TwoColor is a logical argument indicating whether the color map is one
    %or two colors

    if TwoColor 
        nColorMap = size(ColorMap,1);
        MiddleIdx = round(nColorMap/2);

        if mod(nhues,2)==0
            C = zeros(nhues+2,3);

            for i = 1:3
                C(1:round(nhues/2)+1,i) = linspace(ColorMap(1,i),ColorMap(MiddleIdx,i),round(nhues/2)+1)';
                C(round(nhues/2)+2:end,i) = linspace(ColorMap(MiddleIdx,i),ColorMap(end,i),round(nhues/2)+1)';
            end

            C(round(nhues/2)+1:round(nhues/2)+2,:) = [];
        else
            C = zeros(nhues+1,3);

            for i = 1:3
                C(1:round(nhues/2),i) = linspace(ColorMap(1,i),ColorMap(MiddleIdx,i),round(nhues/2)+1)';
                C(round(nhues/2)+1:end,i) = linspace(ColorMap(MiddleIdx,i),ColorMap(end,i),round(nhues/2)+1)';
            end

            C(round(nhues/2)+1,:) = [];
        end
    else
        C = zeros(nhues,3);
        
        for i = 1:3
            C(:,i) = linspace(ColorMap(1,i),ColorMap(end,i),nhues);
        end
    end
end