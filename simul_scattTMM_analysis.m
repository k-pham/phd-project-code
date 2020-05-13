%% parameters

file_dir_data = 'D:\PROJECT\data\simulations\scattTMM\random with water hole 40 80 - fine step size\';
file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\simulations\scattTMM\fine step size\';

scattering_type = 'random';      % options: 'random', 'points'

c_ranges   = 0:10:150;
rho_ranges = 0:10:100;

c0 = 1500;
rho0 = 1000;
c_hole     = 1500;
rho_hole   = 1000;

% set up arrays for saving results of COMPLETE  C & RHO COMPARISON
% scatter_hole_mean_ar = zeros(length(rho_ranges),length(c_ranges));
% scatter_hole_std_ar  = zeros(length(rho_ranges),length(c_ranges));
% scatter_stmm_mean_ar = zeros(length(rho_ranges),length(c_ranges));
% scatter_stmm_std_ar  = zeros(length(rho_ranges),length(c_ranges));
% scatSNR_ar           = zeros(length(rho_ranges),length(c_ranges));
% scatCNR_ar           = zeros(length(rho_ranges),length(c_ranges));


%% loop starts for c & rho
for idx_c = 5 %1:length(c_ranges)
    for idx_r = 9 %1:length(rho_ranges)
        
        % close all
        
        c_scatt = c_ranges(idx_c);
        rho_scatt = rho_ranges(idx_r);
        disp(['STMM ' scattering_type ' c ' num2str(c_scatt) ' rho ' num2str(rho_scatt)])
            %pause
        
        
        %% load data from simu, sensor and recon
        
        file_name = [scattering_type '_SCATT_c' num2str(c_scatt) '_rho' num2str(rho_scatt) ...
                                     '_HOLE_c' num2str(c_hole) '_rho' num2str(rho_hole) ];
        
        SNS = load([file_dir_data file_name '_sensor_data.mat']);
        IMG = load([file_dir_data file_name '_image_data.mat']);
        
        
        %% assign loaded data to structures
        
        sensor.data   = SNS.sensor_data;
        sensor.params = SNS.params;
        simu          = SNS.simu;
        
        sensor.kgrid    = kWaveGrid(sensor.params.Nx, sensor.params.dx, sensor.params.Ny, sensor.params.dy);
        sensor.kgrid.dt = simu.kgrid.dt;
        sensor.kgrid.Nt = simu.kgrid.Nt;
        sensor.t_array  = simu.kgrid.t_array;
        
        image.data    = squeeze(IMG.volume_data);
        image.spacing = IMG.volume_spacing;
        image.kgrid   = IMG.kgrid;
        image.t_array = IMG.t_array;
        image.c0      = c0;
        
        clear SNS IMG
        
        
        %% look at simulation masks
        
        % mask = simu.medium.density;
        % 
        % figure
        % imagesc(simu.kgrid.x_vec*1e3,simu.kgrid.y_vec*1e3,mask')
        %     axis image
        %     xlabel('x position / mm')
        %     ylabel('y position / mm')
        %     colorbar
        
        
        %% look_at_central_sensor_data
        
        % x_centre = round(sensor.params.Nx/2);
        % 
        % figure
        % plot(sensor.t_array*1e6,sensor.data(x_centre,:)) % (50:end)
        %     title('unfiltered central sensor data')
        %     xlabel('time / \mus')
        %     ylabel('acoustic pressure / Pa')
        % 
        % figure
        % imagesc(sensor.kgrid.x_vec*1e3,sensor.t_array(50:end)*1e6,sensor.data(:,50:end)')    
        
        
        %% receive beam forming to look at A-line data
        % average consecutive sensor points to make larger effective sensor
        % element, with more directional response (i.e. A-line)
        
        % element_size   = 11;          % number of sensor elements to average
        % element_centre = round(sensor.params.Nx/2);
        % element_bounds = round(element_centre-element_size/2+0.5 : element_centre+element_size/2-0.5);
        % 
        % element_data   = sensor.data(element_bounds,:);
        % 
        % % element_weights = ones(element_size,1);     % constant weights
        % element_weights = getWin(element_size,'Gaussian');      % gaussian window weights
        % element_weights = element_weights / sum(element_weights);
        % 
        % element_bmform = element_weights .* element_data;
        % element_bmform = sum(element_bmform,1);
        % 
        % figure
        % plot(sensor.t_array*1e6,element_bmform)
        %     title([num2str(element_size) ' sensor points averaged (Gaussian weights): ' num2str(element_bounds(1)) ' - ' num2str(element_bounds(end))])
        %     xlabel('time / \mus')
        %     ylabel('acoustic pressure / Pa')
        
        
        %% frequency filter data before recon
        
        centre_freqs = (1:1:35)*1e6;
        bandwidths   = (1:1:40)*1e6;
        
        % set up arrays for saving results of COMPLETE F & BW TEST
        scatter_hole_mean_ar = zeros(length(bandwidths),length(centre_freqs));
        scatter_hole_std_ar  = zeros(length(bandwidths),length(centre_freqs));
        scatter_stmm_mean_ar = zeros(length(bandwidths),length(centre_freqs));
        scatter_stmm_std_ar  = zeros(length(bandwidths),length(centre_freqs));
        scatSNR_ar           = zeros(length(bandwidths),length(centre_freqs));
        scatCNR_ar           = zeros(length(bandwidths),length(centre_freqs));
        
        for cf_index = 1:length(centre_freqs)
            for bw_index = 1:length(bandwidths)

                centre_freq  = centre_freqs(cf_index);
                bandwidth    = bandwidths(bw_index);
                bandwidth_pc = bandwidth / centre_freq * 100;
                disp(['FILTER ' num2str(centre_freq/1e6) ' bw ' num2str(bandwidth/1e6)])
                    %pause
                
                % apply frequency filter (smooth edges of time series first)
                win = getWin(sensor.kgrid.Nt,'Tukey');
                sensor_smoothedge.data = sensor.data .* win';
                sensor_filtered.data = gaussianFilter(sensor_smoothedge.data,1/sensor.kgrid.dt,centre_freq,bandwidth_pc,false);
                
                % plot central sensor data unfiltered/smoothedged/filtered
                x_centre_sensor = round(sensor.params.Nx/2);
                x_centre_source = round(size(simu.source.p,1)/2);
                figure(1)
                clf(1)
                hold on
                plot(sensor.t_array*1e6,sensor.data(x_centre_sensor,:),'b')
                plot(sensor.t_array*1e6,sensor_smoothedge.data(x_centre_sensor,:),'g')
                plot(sensor.t_array*1e6,sensor_filtered.data(x_centre_sensor,:),'r')
                    title([scattering_type ' c ' num2str(c_scatt) ' rho ' num2str(rho_scatt) ',' ...
                                           ' freq ' num2str(centre_freq/1e6) ' bw ' num2str(bandwidth/1e6)])
                    xlabel('time / \mus')
                    ylabel('acoustic pressure / Pa')
                    legend('sensor data - unfiltered','sensor data - smooth edge','sensor data - filtered')
                    axis([0,10.24,-0.5,0.5])
                
                % evaluate frequency spectra of unfiltered/filtered sensor data & source
                [sensor.freq_axis,            sensor.freq_data            ] = spect(sensor.data(x_centre_sensor,:),            1/sensor.kgrid.dt);
                [sensor_smoothedge.freq_axis, sensor_smoothedge.freq_data ] = spect(sensor_smoothedge.data(x_centre_sensor,:), 1/sensor.kgrid.dt);
                [sensor_filtered.freq_axis,   sensor_filtered.freq_data   ] = spect(sensor_filtered.data(x_centre_sensor,:),   1/sensor.kgrid.dt);
                [simu.source.freq_axis,       simu.source.freq_data       ] = spect(simu.source.p(x_centre_source,:),          1/simu.kgrid.dt  );
                
                % plot frequency spectra of unfiltered/filtered sensor data & source
                figure(2)
                clf(2)
                hold on
                plot(sensor.freq_axis/1e6,            sensor.freq_data,            'b')
                plot(sensor_smoothedge.freq_axis/1e6, sensor_smoothedge.freq_data, 'g')
                plot(sensor_filtered.freq_axis/1e6,   sensor_filtered.freq_data,   'r')
                plot(simu.source.freq_axis/1e6,       simu.source.freq_data,       'k')
                    title([scattering_type ' c ' num2str(c_scatt) ' rho ' num2str(rho_scatt) ',' ...
                                           ' freq ' num2str(centre_freq/1e6) ' bw ' num2str(bandwidth/1e6)])
                    xlim([0,100])
                    xlabel('frequency / MHz')
                    ylabel('amplitude')
                    legend('sensor data - unfiltered','sensor data - smooth edge','sensor data - filtered','source')
                
                % recon with narrowband data and show image
                image_filtered.data = reconstruct2dUSimage(sensor_filtered.data, sensor.params, c0);
                
                % plot filtered recon image
                fig_img = figure(3);
                clf(3)
                imagesc(image.kgrid.x_vec*1e3,image.t_array*c0*1e3,image_filtered.data')
                    axis image
                    title([scattering_type ' c ' num2str(c_scatt) ' rho ' num2str(rho_scatt) ',' ...
                                           ' freq ' num2str(centre_freq/1e6) ' bw ' num2str(bandwidth/1e6)])
                    xlabel('x position / mm')
                    ylabel('y position / mm')
                    colorbar
                
                % convert filtered image to struct for quality metric assessment
                image_filtered.spacing = image.spacing;
                image_filtered.kgrid   = image.kgrid;
                image_filtered.t_array = image.t_array;
                image_filtered.c0      = c0;
                
                
                %% plot image

                % figure
                % imagesc(image.kgrid.x_vec*1e3,image.t_array*c0*1e3,image.data')
                %     axis image
                %     title([scattering_type ' c ' num2str(c_scatt) ' rho ' num2str(rho_scatt)])
                %     xlabel('x position / mm')
                %     ylabel('y position / mm')
                %     colorbar
                
                
                %% scattering distributions in hole & stmm, plot if wanted
                
                plot_toggle = true;
                
                if plot_toggle == true
                    fig_distr = figure(4);
                        clf(4)
                        title('scattering distributions')
                        hold on
                        xlabel('pixel intensity')
                        ylabel('count')
                end
                
                [scatter_hole_mean, scatter_hole_std] = get_scattering_distr_in_hole(image_filtered, plot_toggle);
                [scatter_stmm_mean, scatter_stmm_std] = get_scattering_distr_in_tmm(image_filtered, plot_toggle);
                
                if plot_toggle == true
                    legend(gca,'show')
                end
                
                
                %% image quality metrics
                
                scatSNR = scatter_stmm_mean / scatter_hole_mean;
                scatCNR = (scatter_stmm_mean - scatter_hole_mean) / (scatter_stmm_std + scatter_hole_std);
                
                disp('  scatSNR   scatCNR')
                disp([scatSNR,scatCNR])
                
                
                %% put scattering distribution properties and image qual metrics in array for COMPLETE F & BW TEST
                
                scatter_hole_mean_ar(bw_index,cf_index) = scatter_hole_mean;
                scatter_hole_std_ar(bw_index,cf_index)  = scatter_hole_std;
                scatter_stmm_mean_ar(bw_index,cf_index) = scatter_stmm_mean;
                scatter_stmm_std_ar(bw_index,cf_index)  = scatter_stmm_std;
                
                scatSNR_ar(bw_index,cf_index) = scatSNR;
                scatCNR_ar(bw_index,cf_index) = scatCNR;
                
                
                %% save images for COMPLETE F & BW TEST
                
                saveas(fig_img  ,[file_dir_figs 'freq filtering complete f bw test/f' num2str(centre_freq/1e6) '_bw' num2str(bandwidth/1e6) '_image.jpg'])
                saveas(fig_distr,[file_dir_figs 'freq filtering complete f bw test/f' num2str(centre_freq/1e6) '_bw' num2str(bandwidth/1e6) '_distr.jpg'])
                
                
            end
        end
        
        %% save results of COMPLETE F & BW TEST
        
        save([file_dir_figs 'freq filtering complete f bw test/random_c' num2str(c_scatt) '_rho ' num2str(rho_scatt) '_completeFBWtest.mat'], ...
            'scatter_hole_mean_ar','scatter_hole_std_ar','scatter_stmm_mean_ar','scatter_stmm_std_ar', ...
            'scatSNR_ar','scatCNR_ar')
        
        
        %% put scattering distribution properties and image qual metrics in array for COMPLETE  C & RHO COMPARISON
        
        % scatter_hole_mean_ar(idx_r,idx_c) = scatter_hole_mean;
        % scatter_hole_std_ar(idx_r,idx_c)  = scatter_hole_std;
        % scatter_stmm_mean_ar(idx_r,idx_c) = scatter_stmm_mean;
        % scatter_stmm_std_ar(idx_r,idx_c)  = scatter_stmm_std;
        % 
        % scatSNR_ar(idx_r,idx_c) = scatSNR;
        % scatCNR_ar(idx_r,idx_c) = scatCNR;
        
        
        %% save images for COMPLETE  C & RHO COMPARISON
        
        % saveas(fig_img  ,[file_dir_figs 'complete c rho scattering comparison/c' num2str(c_scatt) '_rho' num2str(rho_scatt) '_image.jpg'])
        % saveas(fig_distr,[file_dir_figs 'complete c rho scattering comparison/c' num2str(c_scatt) '_rho' num2str(rho_scatt) '_distr.jpg'])
        
        
%% loop ends for c & rho
    end
end

%% save results of COMPLETE  C & RHO COMPARISON

% save([file_dir_figs 'complete c rho scattering comparison/completeCRHOtest_f' num2str(centre_freq/1e6) '_bw' num2str(bandwidth/1e6) '.mat'], ...
%     'scatter_hole_mean_ar','scatter_hole_std_ar','scatter_stmm_mean_ar','scatter_stmm_std_ar', ...
%     'scatSNR_ar','scatCNR_ar')


%% plot results of COMPLETE F & BW TEST

centre_freqs = 1:1:35;
bandwidths   = 1:1:40;
figure('Position',[ 100,550,560,420]), imagesc(centre_freqs,bandwidths,scatter_hole_mean_ar), title('scatter hole mean'), xlabel('centre freq / MHz'), ylabel('bandwidth / MHz'), colorbar
figure('Position',[ 100, 30,560,420]), imagesc(centre_freqs,bandwidths,scatter_hole_std_ar) , title('scatter hole std') , xlabel('centre freq / MHz'), ylabel('bandwidth / MHz'), colorbar
figure('Position',[ 700,550,560,420]), imagesc(centre_freqs,bandwidths,scatter_stmm_mean_ar), title('scatter stmm mean'), xlabel('centre freq / MHz'), ylabel('bandwidth / MHz'), colorbar
figure('Position',[ 700, 30,560,420]), imagesc(centre_freqs,bandwidths,scatter_stmm_std_ar) , title('scatter stmm std') , xlabel('centre freq / MHz'), ylabel('bandwidth / MHz'), colorbar
figure('Position',[1300,550,560,420]), imagesc(centre_freqs,bandwidths,scatSNR_ar)          , title('scatSNR')          , xlabel('centre freq / MHz'), ylabel('bandwidth / MHz'), colorbar, caxis([1,3.3])
figure('Position',[1300, 30,560,420]), imagesc(centre_freqs,bandwidths,scatCNR_ar)          , title('scatCNR')          , xlabel('centre freq / MHz'), ylabel('bandwidth / MHz'), colorbar, caxis([0,1.05])


%% plot results of COMPLETE C & RHO COMPARISON

c_ranges   = 0:10:150;
rho_ranges = 0:10:100;
figure('Position',[ 100,550,560,420]), imagesc(c_ranges,rho_ranges,scatter_hole_mean_ar), title('scatter hole mean'), xlabel('c range / m/s'), ylabel('\rho range / kg/m3'), colorbar
figure('Position',[ 100, 30,560,420]), imagesc(c_ranges,rho_ranges,scatter_hole_std_ar) , title('scatter hole std') , xlabel('c range / m/s'), ylabel('\rho range / kg/m3'), colorbar
figure('Position',[ 700,550,560,420]), imagesc(c_ranges,rho_ranges,scatter_stmm_mean_ar), title('scatter stmm mean'), xlabel('c range / m/s'), ylabel('\rho range / kg/m3'), colorbar
figure('Position',[ 700, 30,560,420]), imagesc(c_ranges,rho_ranges,scatter_stmm_std_ar) , title('scatter stmm std') , xlabel('c range / m/s'), ylabel('\rho range / kg/m3'), colorbar
figure('Position',[1300,550,560,420]), imagesc(c_ranges,rho_ranges,scatSNR_ar)          , title('scatSNR')          , xlabel('c range / m/s'), ylabel('\rho range / kg/m3'), colorbar
figure('Position',[1300, 30,560,420]), imagesc(c_ranges,rho_ranges,scatCNR_ar)          , title('scatCNR')          , xlabel('c range / m/s'), ylabel('\rho range / kg/m3'), colorbar


%% compare with experimental data
% 181204 atmm with orgasol

file_dir = '..\data\imagingUS\';
file_name = '181204\atmm_orgasol1_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_03s08m21h_04-12-18_avg1_2D_raw.SGL';
trigger_delay = 0;
samples_cut_off = 50;
samples_t0_correct = -6;
c0 = 1544;

[exp.data, exp.params] = loadSGL([file_dir file_name]);
exp.kgrid   = kWaveGrid(exp.params.Nx, exp.params.dx, exp.params.Ny ,exp.params.dy);
exp.t_array = linspace(1, exp.params.Nt, exp.params.Nt) * exp.params.dt;

figure
imagesc(exp.kgrid.x_vec*1e3, exp.t_array*1e6, squeeze(exp.data(:,round(exp.params.Nx/2),:))')
    xlabel('x / mm')
    ylabel('time / \mus')
    ylim([0,10])
figure
plot(exp.t_array*1e6, squeeze(exp.data(round(exp.params.Nx/2),round(exp.params.Ny/2),:)))
    xlabel('time / \mus')
    ylabel('amplitude')
    xlim([0,10])
    ylim([-0.3,0.2])

figure
hold on
[frequency, f_series] = spect(squeeze(exp.data(75,75,:)),1/exp.params.dt);
plot(frequency/1e6,f_series,'b')
[frequency, f_series] = spect(squeeze(exp.data(75,75,:)),1/exp.params.dt,'Window','Tukey'); 
plot(frequency/1e6,f_series,'g')
    xlim([0,60])
    ylim([0,0.01])
    xlabel('frequency / MHz')
    ylabel('amplitude')


%% check noise statistics of experimental data
% 191126 resolution27umPlanar in 2D // 181119 hair in 2D // 

file_dir = '..\data\imagingUS\';
% file_name = '191126\resolution27umPlanar_BK31[CNT]_trolley_scrambled_2D@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_13s20m12h_26-11-19_avg1_2D_raw.SGL';
% file_name = '181119\hair_BK31[CNT]_6@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_17s24m17h_19-11-18_avg1_2D_raw.SGL';
file_name = '180625\polymerLeafFlat_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_08s15m19h_25-06-18_avg1_2D_raw.SGL';

[exp.data, exp.params] = loadSGL([file_dir file_name]);

% ROI = exp.data(45,60:90,50:500);
ROI = exp.data(70,1:40,1:180);
% ROI = exp.data(70,:,1:1000);

figure, imagesc(squeeze(ROI))

binwidth   = max(ROI(:))/100;
binedges   = min(ROI(:)) : binwidth : max(ROI(:));
bincount   = histcounts(ROI, binedges);
bincentres = 0.5 * ( binedges(2:end) + binedges(1:end-1) );

figure, plot(bincentres,bincount)
hold on

fittttt = fit(bincentres', bincount', 'gauss1');
gcf, plot(fittttt)


%% LOCAL FUNCTIONS

function ROI = get_discROI_at_hole_in_image(image, fractional_radius)

    hole_x = 0;
    hole_y = 2.35e-3;
    hole_radius = 0.5e-3;
    
    vec_x = image.kgrid.x_vec;
    vec_y = image.t_array*image.c0;
    
    distance = sqrt((vec_x - hole_x).^2 + (vec_y - hole_y).^2);
    
    ROI = distance < hole_radius * fractional_radius;
    
end

function [scatter_hole_mean, scatter_hole_std] = get_scattering_distr_in_hole(image, plot_toggle)

    mask = get_discROI_at_hole_in_image(image, 0.9);
    
    ROI = image.data(mask);
    
    scatter_hole_mean = mean(ROI(:));
    scatter_hole_std  = std(ROI(:));
    
    if plot_toggle == true
        plot_histogram_of_scattering_distr(ROI, 'scatter hole', 'Normalise', true)
    end
    
end

function [scatter_tmm_mean, scatter_tmm_std] = get_scattering_distr_in_tmm(image, plot_toggle)

    hole_notoutside = get_discROI_at_hole_in_image(image, 1.1);
    hole_outside    = not(hole_notoutside);
    large_hole      = get_discROI_at_hole_in_image(image, 2);
    
    mask = and(hole_outside,large_hole);
    
    ROI = image.data(mask);
    
    scatter_tmm_mean = mean(ROI(:));
    scatter_tmm_std  = std(ROI(:));
    
    if plot_toggle == true
        plot_histogram_of_scattering_distr(ROI, 'scatter tmm', 'Normalise', true)
    end
    
end

function plot_histogram_of_scattering_distr(ROI, legend_entry, varargin)

    % set default
    num_req_input_variables = 2;
    toNormalise = false;
    
    if nargin < num_req_input_variables
        error('Incorrect number of inputs.');
    elseif ~isempty(varargin)
        for input_index = 1:2:length(varargin)
            switch varargin{input_index}
                case 'Normalise'
                    toNormalise = varargin{input_index + 1};
                otherwise
                    error('Unknown optional input.');
            end
        end
    end
    
    binwidth   = max(ROI(:))/100;    
    binedges   = 0:binwidth:max(ROI(:))*1.1+binwidth;
    bincount   = histcounts(ROI,binedges);
    bincentres = 0.5*(binedges(2:end)+binedges(1:end-1));
	
    if toNormalise
        bincount = bincount / max(bincount);
    end
    
    figure(gcf)
    % histogram(ROI,'BinWidth',10,'DisplayStyle','stairs','DisplayName',legend_entry)
    plot(bincentres,bincount,'DisplayName',legend_entry)
    
end


