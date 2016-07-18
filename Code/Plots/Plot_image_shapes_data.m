%% Plot out of bag error vs n on mnist dataset

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load image_shapes_data

X0 = X_image(:,:,Y==0);
X1 = X_image(:,:,Y==1);
Y0 = Y(Y==0);
Y1 = Y(Y==1);

for i = 1:3
    subplot(2,3,i)
    imshow(X0(:,:,i))
    if i == 1
        ylabel('Class 0','FontSize',14)
    end
    subplot(2,3,i+3)
    imshow(X1(:,:,i))
    if i == 1
        ylabel('Class 1','FontSize',14)
    end
end
save_fig(gcf,[rerfPath 'RandomerForest/Figures/Image_shapes_data'])