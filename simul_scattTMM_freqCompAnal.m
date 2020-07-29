% load simu
load('D:\PROJECT\data\simulations\scattTMM\freq compounding & coherent filter test2\random_SCATT_c40_rho80_hole_OBJECT_c1500_rho1000_x768_y256_simu.mat')

% load various images
load('D:\PROJECT\data\simulations\scattTMM\freq compounding & coherent filter test2\random_SCATT_c40_rho80_hole_OBJECT_c1500_rho1000_x768_y256_image.mat')
image_unfiltered = image;

load('D:\PROJECT\data\simulations\scattTMM\freq compounding & coherent filter test2\random_SCATT_c40_rho80_hole_OBJECT_c1500_rho1000_x768_y256_FILTER_f1_bw20_image.mat')
image_gaussian   = image;

load('D:\PROJECT\data\simulations\scattTMM\freq compounding & coherent filter test2\random_SCATT_c40_rho80_hole_OBJECT_c1500_rho1000_x768_y256_COMPOUND_incoherent_f1-10_bw2_image.mat')
image_incoherent = image;

load('D:\PROJECT\data\simulations\scattTMM\freq compounding & coherent filter test2\random_SCATT_c40_rho80_hole_OBJECT_c1500_rho1000_x768_y256_COMPOUND_coherent_f1-10_bw2_image.mat')
image_coherent   = image;

load('D:\PROJECT\data\simulations\scattTMM\freq compounding & coherent filter test2\random_SCATT_c40_rho80_hole_OBJECT_c1500_rho1000_x768_y256_CFILTER_f1-10_bw2_image.mat')
image_cfilter    = image;

clear image

% subtract to get difference images
image_diff1 = image_unfiltered;
image_diff2 = image_unfiltered;

image_diff1.data = image_coherent.data - image_incoherent.data;
image_diff2.data = image_coherent.data - image_cfilter.data;

% plot difference images
fig_imagdiff1 = plot_image_data(image_diff1, simu);
fig_imagdiff2 = plot_image_data(image_diff2, simu);


%% FUNCTIONS

function fig_imag = plot_image_data(image, simu)
% plots:    image.data
% requires: image.data
%           image.kgrid   - for axes
%           image.t_array - for axes
%           simu.params   - for c0, titles

    fig_imag = figure;
    imagesc(image.kgrid.x_vec*1e3, image.t_array*simu.params.c0*1e3, image.data')
                     % omit factor 1/2 in depth because of doubled depth bug
        axis image
        title([simu.params.scatt_type   ' c ' num2str(simu.params.scatt_c)  ' rho ' num2str(simu.params.scatt_rho) ', ' ...
               simu.params.object_shape ' c ' num2str(simu.params.object_c) ' rho ' num2str(simu.params.object_rho) ])
        xlabel('x position / mm')
        ylabel('y position / mm')
        colorbar
    
end