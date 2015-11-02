function X = hyperspiral(r_inc,phi_inc,npoints)
    t = 0:1/npoints:1;
    d = length(phi_inc) + 1;
    X = zeros(npoints,d);
    phi = zeros(npoints,d-1);
    r = transpose(0:r_inc:(npoints-1)*r_inc);
    %r = ones(npoints,1);
    for i = 1:d-1
        phi(:,i) = transpose(0:phi_inc(i):(npoints-1)*phi_inc(i));
    end
    X(:,1) = r.*sin(phi(:,1));
    X(:,2) = r.*cos(phi(:,1));
    if d > 2
        for i = 2:d-1
            X(:,1) = X(:,1).*sin(phi(:,i));
            X(:,2) = X(:,2).*sin(phi(:,i));
        end
        for i = 3:d
            X(:,i) = X(:,i-1).*cos(phi(:,i-1))./cos(phi(:,i-2))./sin(phi(:,i-1));
        end
    end