%Computes and plots bayes error for Trunk as a function of d

clear
close all
clc

dims = round(logspace(log10(2),3,5));
bayes_error = zeros(size(dims));
for i = 1:length(dims)
    d = dims(i);
    d_idx = 1:d;
    mu1 = 1./d_idx;
    mu0 = -1*mu1;
    Sigma = speye(d);
    delta = mu1 - mu0;
    bayes_error(i) = 1 - normcdf(sqrt(delta*Sigma*delta')/2);
end
plot(dims,bayes_error)
xlabel('# ambient dimensions')
ylabel('bayes error')
set(gca,'XScale','log')
%save_fig(gcf,'Trunk/Trunk_bayes_error_vs_d')