clear
close all
clc

x = 0:pi/10000:4*pi;

h = figure;
subplot(2,2,1)
plot(x,sin(x))
axis square
subplot(2,2,2)
plot(x,cos(x))
axis square
subplot(2,2,3)
plot(x,sin(x).^2)
axis square
subplot(2,2,4)
plot(x,cos(x).^2)
axis square

h.PaperUnits = 'normalized';
h.PaperPosition = [0 0 1 1];
h.Units = 'normalized';
h.Position = [0 0 1 1];

saveas(gcf,'posteriors','pdf')