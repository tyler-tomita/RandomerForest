close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')
Colors.rf = ColorMap(1,:);
Colors.rerf = ColorMap(8,:);
Colors.rr_rf = ColorMap(4,:);
LineWidth = 2;
FontSize = .15;
axWidth = 1.5;
axHeight = 1.5;
axLeft = repmat([FontSize*5,FontSize*10+axWidth,FontSize*15+2*axWidth],1,3);
axBottom = [(FontSize*9+2*axHeight)*ones(1,3),(FontSize*6+axHeight)*ones(1,3),FontSize*3*ones(1,3)];
legWidth = axWidth;
legHeight = axHeight;
legLeft = axLeft(end) + axWidth + FontSize/2;
legBottom = axBottom(6);
figWidth = legLeft(end) + legWidth + FontSize*2;
figHeight = axBottom(1) + axHeight + FontSize*2;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

ns = [1,2,5,10,20,50,100];
L_rf = zeros(length(ns),1);
L_rerf = zeros(length(ns),1);

for i = 1:length(ns)
    n = ns(i);
    fprintf('n = %d\n',n)
    a = 0.01;
    vx = 1;
    e = a/100/sqrt(vx);
    mx1_0 = -a;
    mx2_0 = -e;
    mx_0 = [mx1_0,mx2_0];
    mx1_1 = a;
    mx2_1 = e;
    mx_1 = [mx1_1,mx2_1];
    m1 = 1/2*(mx1_0 + mx1_1);
    m2 = 1/2*(mx2_0 + mx2_1);
    v = vx/n;
    funx_0 = @(x1,x2) 1/2/pi/vx.*exp(-1/2.*((x1-mx1_0).^2/vx + (x2-mx2_0).^2/vx));
    funx_1 = @(x1,x2) 1/2/pi/vx.*exp(-1/2.*((x1-mx1_1).^2/vx + (x2-mx2_1).^2/vx));
    funmu_0 = @(mu1,mu2) integral2(funx_0,mu1,mx1_0+10*sqrt(vx),@(x1) -x1+mu1+mu2,mu2)*1/2/pi/v.*exp(-1/2.*((mu1-m1).^2/v + (mu2-m2).^2/v));
    funmu_0v = @(mu1,mu2) arrayfun(funmu_0,mu1,mu2);
    funmu_1 = @(mu1,mu2) integral2(funx_1,mu1,mx1_1+10*sqrt(vx),@(x1) -x1+mu1+mu2,mu2)*1/2/pi/v.*exp(-1/2.*((mu1-m1).^2/v + (mu2-m2).^2/v));
    funmu_1v = @(mu1,mu2) arrayfun(funmu_1,mu1,mu2);
    
    L_rf(i) = 1/2*normcdf(mx1_0/sqrt(vx)/sqrt(1+1/n))*normcdf(mx2_0/sqrt(vx)/sqrt(1+1/n)) + ...
        1/2*(1 - normcdf(mx1_1/sqrt(vx)/sqrt(1+1/n))*normcdf(mx2_1/sqrt(vx)/sqrt(1+1/n)));
        
    L_rerf(i) = 1/2*(normcdf(mx1_0/sqrt(vx)/sqrt(1+1/n))*normcdf(mx2_0/sqrt(vx)/sqrt(1+1/n)) + ...
        integral2(funmu_0v,-4*sqrt(v),4*sqrt(v),-4*sqrt(v),4*sqrt(v))) + ...
        1/2*(1 - (normcdf(mx1_1/sqrt(vx)/sqrt(1+1/n))*normcdf(mx2_1/sqrt(vx)/sqrt(1+1/n)) + ...
        integral2(funmu_1v,-4*sqrt(v),4*sqrt(v),-4*sqrt(v),4*sqrt(v))));
end

[X1,X2] = meshgrid(-5:.01:5,-5:.01:5);
f_X = 1/2*(1/2/pi/vx.*exp(-1/2.*((X1-mx1_0).^2/vx + (X2-mx2_0).^2/vx)) + 1/2/pi/vx.*exp(-1/2.*((X1-mx1_1).^2/vx + (X2-mx2_1).^2/vx)));
j = 1;
ax(j) = axes;
mesh(X1,X2,f_X)
grid off
xlabel('x_1')
ylabel('x_2')
zlabel('f_X(x)')
ax(j).LineWidth = LineWidth;
ax(j).FontUnits = 'inches';
ax(j).FontSize = FontSize;
ax(j).Units = 'inches';
ax(j).Position = [axLeft(j) axBottom(j) axWidth axHeight];

j = 2;
ax(j) = axes;
heat_map(X1,X2,f_X,parula);
xlabel('x_1')
ylabel('x_2')
ax(j).LineWidth = LineWidth;
ax(j).FontUnits = 'inches';
ax(j).FontSize = FontSize;
ax(j).Units = 'inches';
ax(j).Position = [axLeft(j) axBottom(j) axWidth axHeight];

D = [mx1_1-mx1_0;mx2_1-mx2_0];
S = diag([vx,vx]);
L_bayes = normcdf(-1/2*sqrt(D'*S*D));
j = 3;
ax(j) = axes;
plot(ns,L_rf,ns,L_rerf,[ns(1) ns(end)],L_bayes*ones(1,2),'k','LineWidth',2)
xlabel('n')
ylabel('E[L_n]')
ax(j).LineWidth = LineWidth;
ax(j).FontUnits = 'inches';
ax(j).FontSize = FontSize;
ax(j).Units = 'inches';
ax(j).Position = [axLeft(j) axBottom(j) axWidth axHeight];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ns = [1,2,5,10,20,50,100];
L_rf = zeros(length(ns),1);
L_rerf = zeros(length(ns),1);

for i = 1:length(ns)
    n = ns(i);
    fprintf('n = %d\n',n)
    mx_0r = rotate2d(mx_0,pi/4);
    mx1_0 = mx_0r(1);
    mx2_0 = mx_0r(2);
    mx_1r = rotate2d(mx_1,pi/4);
    mx1_1 = mx_1r(1);
    mx2_1 = mx_1r(2);
    m1 = 1/2*(mx1_0 + mx1_1);
    m2 = 1/2*(mx2_0 + mx2_1);
    v = vx/n;
    funx_0 = @(x1,x2) 1/2/pi/vx.*exp(-1/2.*((x1-mx1_0).^2/vx + (x2-mx2_0).^2/vx));
    funx_1 = @(x1,x2) 1/2/pi/vx.*exp(-1/2.*((x1-mx1_1).^2/vx + (x2-mx2_1).^2/vx));
    funmu_0 = @(mu1,mu2) integral2(funx_0,mu1,mx1_0+10*sqrt(vx),@(x1) -x1+mu1+mu2,mu2)*1/2/pi/v.*exp(-1/2.*((mu1-m1).^2/v + (mu2-m2).^2/v));
    funmu_0v = @(mu1,mu2) arrayfun(funmu_0,mu1,mu2);
    funmu_1 = @(mu1,mu2) integral2(funx_1,mu1,mx1_1+10*sqrt(vx),@(x1) -x1+mu1+mu2,mu2)*1/2/pi/v.*exp(-1/2.*((mu1-m1).^2/v + (mu2-m2).^2/v));
    funmu_1v = @(mu1,mu2) arrayfun(funmu_1,mu1,mu2);
    
    L_rf(i) = 1/2*normcdf(mx1_0/sqrt(vx)/sqrt(1+1/n))*normcdf(mx2_0/sqrt(vx)/sqrt(1+1/n)) + ...
        1/2*(1 - normcdf(mx1_1/sqrt(vx)/sqrt(1+1/n))*normcdf(mx2_1/sqrt(vx)/sqrt(1+1/n)));
        
    L_rerf(i) = 1/2*(normcdf(mx1_0/sqrt(vx)/sqrt(1+1/n))*normcdf(mx2_0/sqrt(vx)/sqrt(1+1/n)) + ...
        integral2(funmu_0v,-4*sqrt(v),4*sqrt(v),-4*sqrt(v),4*sqrt(v))) + ...
        1/2*(1 - (normcdf(mx1_1/sqrt(vx)/sqrt(1+1/n))*normcdf(mx2_1/sqrt(vx)/sqrt(1+1/n)) + ...
        integral2(funmu_1v,-4*sqrt(v),4*sqrt(v),-4*sqrt(v),4*sqrt(v))));
end

[X1,X2] = meshgrid(-5:.01:5,-5:.01:5);
f_X = 1/2*(1/2/pi/vx.*exp(-1/2.*((X1-mx1_0).^2/vx + (X2-mx2_0).^2/vx)) + 1/2/pi/vx.*exp(-1/2.*((X1-mx1_1).^2/vx + (X2-mx2_1).^2/vx)));
j = 4;
ax(j) = axes;
mesh(X1,X2,f_X)
grid off
xlabel('x_1')
ylabel('x_2')
zlabel('f_X(x)')
ax(j).LineWidth = LineWidth;
ax(j).FontUnits = 'inches';
ax(j).FontSize = FontSize;
ax(j).Units = 'inches';
ax(j).Position = [axLeft(j) axBottom(j) axWidth axHeight];

j = 5;
ax(j) = axes;
heat_map(X1,X2,f_X,parula);
xlabel('x_1')
ylabel('x_2')
ax(j).LineWidth = LineWidth;
ax(j).FontUnits = 'inches';
ax(j).FontSize = FontSize;
ax(j).Units = 'inches';
ax(j).Position = [axLeft(j) axBottom(j) axWidth axHeight];

D = [mx1_1-mx1_0;mx2_1-mx2_0];
S = diag([vx,vx]);
L_bayes = normcdf(-1/2*sqrt(D'*S*D));
j = 6;
ax(j) = axes;
plot(ns,L_rf,ns,L_rerf,[ns(1) ns(end)],L_bayes*ones(1,2),'k','LineWidth',2)
xlabel('n')
ylabel('E[L_n]')
ax(j).LineWidth = LineWidth;
ax(j).FontUnits = 'inches';
ax(j).FontSize = FontSize;
ax(j).Units = 'inches';
ax(j).Position = [axLeft(j) axBottom(j) axWidth axHeight];

[lh,objh] = legend('RF','RerF','Bayes');
lh.Units = 'inches';
lh.Position = [legLeft legBottom legWidth legHeight];
lh.Box = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ns = [1,2,5,10,20,50,100];
L_rf = zeros(length(ns),1);
L_rerf = zeros(length(ns),1);

for i = 1:length(ns)
    n = ns(i);
    fprintf('n = %d\n',n)
    mx_0r = rotate2d(mx_0,pi/8);
    mx1_0 = mx_0r(1);
    mx2_0 = mx_0r(2);
    mx_1r = rotate2d(mx_1,pi/8);
    mx1_1 = mx_1r(1);
    mx2_1 = mx_1r(2);
    m1 = 1/2*(mx1_0 + mx1_1);
    m2 = 1/2*(mx2_0 + mx2_1);
    v = vx/n;
    funx_0 = @(x1,x2) 1/2/pi/vx.*exp(-1/2.*((x1-mx1_0).^2/vx + (x2-mx2_0).^2/vx));
    funx_1 = @(x1,x2) 1/2/pi/vx.*exp(-1/2.*((x1-mx1_1).^2/vx + (x2-mx2_1).^2/vx));
    funmu_0 = @(mu1,mu2) integral2(funx_0,mu1,mx1_0+10*sqrt(vx),@(x1) -x1+mu1+mu2,mu2)*1/2/pi/v.*exp(-1/2.*((mu1-m1).^2/v + (mu2-m2).^2/v));
    funmu_0v = @(mu1,mu2) arrayfun(funmu_0,mu1,mu2);
    funmu_1 = @(mu1,mu2) integral2(funx_1,mu1,mx1_1+10*sqrt(vx),@(x1) -x1+mu1+mu2,mu2)*1/2/pi/v.*exp(-1/2.*((mu1-m1).^2/v + (mu2-m2).^2/v));
    funmu_1v = @(mu1,mu2) arrayfun(funmu_1,mu1,mu2);
    
    L_rf(i) = 1/2*normcdf(mx1_0/sqrt(vx)/sqrt(1+1/n))*normcdf(mx2_0/sqrt(vx)/sqrt(1+1/n)) + ...
        1/2*(1 - normcdf(mx1_1/sqrt(vx)/sqrt(1+1/n))*normcdf(mx2_1/sqrt(vx)/sqrt(1+1/n)));
        
    L_rerf(i) = 1/2*(normcdf(mx1_0/sqrt(vx)/sqrt(1+1/n))*normcdf(mx2_0/sqrt(vx)/sqrt(1+1/n)) + ...
        integral2(funmu_0v,-4*sqrt(v),4*sqrt(v),-4*sqrt(v),4*sqrt(v))) + ...
        1/2*(1 - (normcdf(mx1_1/sqrt(vx)/sqrt(1+1/n))*normcdf(mx2_1/sqrt(vx)/sqrt(1+1/n)) + ...
        integral2(funmu_1v,-4*sqrt(v),4*sqrt(v),-4*sqrt(v),4*sqrt(v))));
end

[X1,X2] = meshgrid(-5:.01:5,-5:.01:5);
f_X = 1/2*(1/2/pi/vx.*exp(-1/2.*((X1-mx1_0).^2/vx + (X2-mx2_0).^2/vx)) + 1/2/pi/vx.*exp(-1/2.*((X1-mx1_1).^2/vx + (X2-mx2_1).^2/vx)));
j = 7;
ax(j) = axes;
mesh(X1,X2,f_X)
grid off
xlabel('x_1')
ylabel('x_2')
zlabel('f_X(x)')
ax(j).LineWidth = LineWidth;
ax(j).FontUnits = 'inches';
ax(j).FontSize = FontSize;
ax(j).Units = 'inches';
ax(j).Position = [axLeft(j) axBottom(j) axWidth axHeight];

j = 8;
ax(j) = axes;
heat_map(X1,X2,f_X,parula);
xlabel('x_1')
ylabel('x_2')
ax(j).LineWidth = LineWidth;
ax(j).FontUnits = 'inches';
ax(j).FontSize = FontSize;
ax(j).Units = 'inches';
ax(j).Position = [axLeft(j) axBottom(j) axWidth axHeight];

D = [mx1_1-mx1_0;mx2_1-mx2_0];
S = diag([vx,vx]);
L_bayes = normcdf(-1/2*sqrt(D'*S*D));
j = 9;
ax(j) = axes;
plot(ns,L_rf,ns,L_rerf,[ns(1) ns(end)],L_bayes*ones(1,2),'k','LineWidth',2)
xlabel('n')
ylabel('E[L_n]')
ax(j).LineWidth = LineWidth;
ax(j).FontUnits = 'inches';
ax(j).FontSize = FontSize;
ax(j).Units = 'inches';
ax(j).Position = [axLeft(j) axBottom(j) axWidth axHeight];

% save_fig(gcf,[rerfPath 'RandomerForest/Figures/2017.03.12/Stump_ensemble_theoretical_error'],{'fig','pdf','png'})