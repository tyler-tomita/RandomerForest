function [X,Y] = spiral_scatter(x1_off,x2_off,r_i,r_f,npoints)
phi1_i = 0;
phi1_f = 6*pi;
phi1_inc = (phi1_f - phi1_i)/npoints;
phi1 = phi1_i:phi1_inc:phi1_f;
r_inc = (r_f - r_i)/npoints;
r = r_i:r_inc:r_f;
x11 = r.*cos(phi1);
x21 = r.*sin(phi1);
x12 = -x11;
x22 = -x21;
x1 = cat(2,x11',x21');
x2 = cat(2,x12',x22');
X = cat(1,x1,x2);
X = X + repmat([x1_off x2_off],size(X,1),1);
Y = cat(1,zeros(size(X,1)/2,1),ones(size(X,1)/2,1));
plot(X(Y==0,1),X(Y==0,2),'b.',X(Y==1,1),X(Y==1,2),'r.')
end