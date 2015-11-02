function s = map2size(x,size_min,size_max,map)

    if strcmp(map,'log')
        x = log(x);
    end

    x = normalize(x);
    s = zeros(length(x),1);

    for i = 1:length(x)
        s(i) = round(x(i)*(size_max-size_min) + size_min);
    end
end
