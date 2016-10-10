function X = rand_hypersphere(n,d,R)
%   X = RAND_HYPERSPHERE(RADIUS) samples n points in a d-dimensional
%   hypersphere of radius R

%     X = zeros(n,d);
%     
    % uniformly randomly sample radii and angles
%     Radii = rand(n,1)*R;
%     Phi = rand(n,d-1)*2*pi;
%     
%     X(:,1) = Radii.*sin(Phi(:,1));
%     X(:,2) = Radii.*cos(Phi(:,1));
%     if d > 2
%         for i = 2:d-1
%             X(:,1) = X(:,1).*sin(Phi(:,i));
%             X(:,2) = X(:,2).*sin(Phi(:,i));
%         end
%         for i = 3:d
%             X(:,i) = X(:,i-1).*cos(Phi(:,i-1))./cos(Phi(:,i-2))./sin(Phi(:,i-1));
%         end
%     end

    Vsphere = pi^(d/2)*R^d/gamma(d/2+1);
    Vcube = (2*R)^d;
    nsample = ceil(1.5*n*Vcube/Vsphere);
    if nsample > 10000
        X = [];
        while size(X,1) < n
            x = rand(n,d)*2*R - R;
            Radii = sqrt(sum(x.^2,2));
            x(Radii>R,:) = [];
            X = [X;x];
        end
    else
        X = rand(nsample,d)*2*R - R;
        Radii = sqrt(sum(X.^2,2));
        X(Radii>R,:) = [];
    end
    X = X(1:n,:);
end