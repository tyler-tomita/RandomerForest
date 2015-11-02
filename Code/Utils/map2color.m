function rgb = map2color(x,map)
    if strcmp(map,'log')
        x = log(x);
    end
    x = normalize(x);
    r = zeros(length(x),1);
    g = zeros(length(x),1);
    b = zeros(length(x),1);

    for i = 1:length(x)
        r(i) = max(x(i)/0.5-1,0);
        g(i) = 1 - abs(x(i)/0.5-1);
        b(i) = max(1-x(i)/0.5,0);
    end
    rgb = [r g b];
end
