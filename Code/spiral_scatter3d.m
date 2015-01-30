function [X,Y,Z] = spiral_scatter3d(x1_off,x2_off,x3_off,r_i,r_f,npoints)
phi1_i = 0;
phi1_f = 10*pi;
phi1_inc = (phi1_f - phi1_i)/npoints;
phi1 = phi1_i:phi1_inc:phi1_f;
phi2_i = 0;
phi2_f = 10*pi;
phi2_inc = (phi2_f - phi2_i)/npoints;
phi2 = phi2_i:phi2_inc:phi2_f;
r_inc = (r_f - r_i)/npoints;
r = r_i:r_inc:r_f;
x11 = r.*cos(phi1);
x21 = r.*sin(phi1).*cos(phi2);
x31 = r.*sin(phi1).*sin(phi2);
x12 = -x11;
x22 = -x21;
x32 = -x31;
x1 = cat(2,x11',x21',x31');
x2 = cat(2,x12',x22',x32');
X = cat(1,x1,x2);
X = X + repmat([x1_off x2_off x3_off],size(X,1),1);
Y = cat(1,zeros(size(X,1)/2,1),ones(size(X,1)/2,1));
plot3(X(Y==0,1),X(Y==0,2),X(Y==0,3),'b.',X(Y==1,1),X(Y==1,2),X(Y==1,3),'r.')
end