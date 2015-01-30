function [X] = spiral_scatter3d_jv(x_init,phi1_inc,phi2_inc,r_inc,npoints)
phi1_i = 0;
phi1_f = 4*pi;
phi1 = phi1_i:phi1_inc:phi1_inc*npoints;
phi2_i = 0;
phi2 = phi2_i:phi2_inc:phi2_inc*npoints;
r = 0:r_inc:r_inc*npoints;
x11 = r.*cos(phi1);
x21 = r.*sin(phi1).*cos(phi2);
x31 = x11*0; %r.*sin(phi1).*sin(phi2).*cos(phi2);
% x12 = -x11;
% x22 = -x21;
% x32 = -x31;
x1 = cat(2,x11',x21',x31');
% x2 = cat(2,x12',x22',x32');
% X = cat(1,x1,x2);
X = x1; % + repmat([x_init],size(X,1),1);
% Y = cat(1,zeros(size(X,1)/2,1),ones(size(X,1)/2,1));
% plot3(X(Y==0,1),X(Y==0,2),X(Y==0,3),'b',X(Y==1,1),X(Y==1,2),X(Y==1,3),'r')
plot3(X(:,1),X(:,2),X(:,3),'b')
end