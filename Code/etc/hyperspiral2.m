function X = hyperspiral2(r_final,phi_final,npoints)
    t = transpose(0:1/npoints:1);
    d = length(phi_final) + 1;
    X = zeros(length(t),d);
    phi = zeros(length(t),d-1);
    r = t*r_final;
    %r = ones(npoints,1);
    for i = 1:d-1
        phi(:,i) = t*phi_final(i);
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