function [X,Y] = Spiral_2class(n,d)
%Generates 2 noisy spirals for binary classification

    fpath = mfilename('fullpath');
    findex = strfind(fpath,'/');
    rootDir=fpath(1:findex(end-1));
    p = genpath(rootDir);
    gits=strfind(p,'.git');
    colons=strfind(p,':');
    for i=0:length(gits)-1
    endGit=find(colons>gits(end-i),1);
    p(colons(endGit-1):colons(endGit)-1)=[];
    end
    addpath(p);

    %phi_final = (d-1)*2*pi*randi(d-1,1,d-1);
    if d == 2
        phi_final = 4*pi;
    else
        %phi_final = (d-1)*2*pi*linspace(1,2,d-1);
        phi_final = 4*pi*linspace(1,2,d-1);
    end
    r_final = 1;
    X0 = zeros(round(n/2),d);
    Y0 = zeros(size(X0,1),1);
    X1 = zeros(round(n/2),d);
    Y1 = ones(size(X1,1),1);

    Mu = zeros(1,d);
    Sigma = 0.0075*ones(1,d);
    for i = 1:length(Y0)
        t0 = rand;
        r = t0*r_final;
        phi = t0*phi_final;
        Mu(1) = r.*sin(phi(1));
        Mu(2) = r.*cos(phi(1));
        if d > 2
            for j = 2:d-1
                Mu(1) = Mu(1).*sin(phi(j));
                Mu(2) = Mu(2).*sin(phi(j));
            end
            for j = 3:d
                Mu(j) = Mu(j-1).*cos(phi(j-1))./cos(phi(j-2))./sin(phi(j-1));
            end
        end
        X0(i,:) = mvnrnd(Mu,Sigma,1);

        t1 = rand;
        r = t1*r_final;
        phi = t1*phi_final;
        Mu(1) = r.*sin(phi(1));
        Mu(2) = r.*cos(phi(1));
        if d > 2
            for j = 2:d-1
                Mu(1) = Mu(1).*sin(phi(j));
                Mu(2) = Mu(2).*sin(phi(j));
            end
            for j = 3:d
                Mu(j) = Mu(j-1).*cos(phi(j-1))./cos(phi(j-2))./sin(phi(j-1));
            end
        end
        X1(i,:) = mvnrnd(Mu,Sigma,1);
    end
    X1 = -X1;
    X = cat(1,X0,X1);
    Y = cat(1,Y0,Y1);
end
