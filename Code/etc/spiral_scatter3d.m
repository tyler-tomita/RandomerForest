function X = spiral_scatter3d(phi1_inc,phi2_inc,r_inc,npoints)
phi1 = 0:phi1_inc:phi1_inc*npoints;
phi2 = 0:phi2_inc:phi2_inc*npoints;
%phi1 = ones(1,length(phi1));
phi2 = ones(1,length(phi2));
%phi2 = pi/2;
r = 0:r_inc:r_inc*npoints;
%r = 1;
%[phi1,phi2] = meshgrid(phi1,phi2);
x11 = r.*cos(phi1).*sin(phi2);
%x21 = r.*sin(phi1).*cos(phi2);
%x31 = r.*sin(phi1).*sin(phi2); %r.*sin(phi1).*sin(phi2).*cos(phi2);
x21 = r.*sin(phi1).*sin(phi2);
x31 = r.*cos(phi2);
% x12 = -x11;
% x22 = -x21;
% x32 = -x31;
x1 = cat(2,x11',x21',x31');
% x2 = cat(2,x12',x22',x32');
% X = cat(1,x1,x2);
X = x1; % + repmat([x_init],size(X,1),1);
%Y = cat(1,zeros(size(X,1)/2,1),ones(size(X,1)/2,1));
%plot3(X(Y==0,1),X(Y==0,2),X(Y==0,3),'b',X(Y==1,1),X(Y==1,2),X(Y==1,3),'r')
%surf(x11,x21,x31)
plot3(X(:,1),X(:,2),X(:,3))
end