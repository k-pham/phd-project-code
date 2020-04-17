%% load images

file_dir_data = 'D:\PROJECT\data\simulations\scattTMM\random with water hole\';
file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\simulations\scattTMM\';

scattering_type = 'random';      % options: 'random', 'points'

c_ranges   = 0:10:150;
rho_ranges = 0:10:100;

c0 = 1500;
rho0 = 1000;
c_hole     = 1500;
rho_hole   = 1000;

for idx_c = 5%1:length(c_ranges)
    for idx_r = 9%1:length(rho_ranges)
        
        % close all
        
        c_scatt = c_ranges(idx_c);
        rho_scatt = rho_ranges(idx_r);
        
        
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
        
        x_centre = round(sensor.params.Nx/2);
        
        figure
        plot(sensor.t_array*1e6,sensor.data(x_centre,:)) % (50:end)
            title('unfiltered central sensor data')
            xlabel('time / \mus')
            ylabel('acoustic pressure / Pa')

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
        
        
        %% narrowband data before recon
        
        centre_freq = 10e6;
        bandwidth = 10e6;
        
        bandwidth_pc = bandwidth / centre_freq * 100;
        
        win = getWin(sensor.kgrid.Nt,'Tukey');
        sensor_data_smoothedge = sensor.data .* win';
        sensor_data_filtered = gaussianFilter(sensor_data_smoothedge,1/sensor.kgrid.dt,centre_freq,bandwidth_pc,false);
        
        % plot filtered sensor data
        x_centre = round(sensor.params.Nx/2);
        figure
        imagesc(sensor_data_smoothedge)
            title('smooth-edged sensor data')
        figure
        imagesc(sensor_data_filtered)
            title(['filtered ' num2str(centre_freq/1e6) ' bw ' num2str(bandwidth/1e6) ' central sensor data'])
        figure
        plot(sensor.t_array*1e6,sensor_data_filtered(x_centre,:)) % (50:end)
            title(['filtered ' num2str(centre_freq/1e6) ' bw ' num2str(bandwidth/1e6) ' central sensor data'])
            xlabel('time / \mus')
            ylabel('acoustic pressure / Pa')
        
        % evaluate & plot frequency spectra of unfiltered/filtered sensor data & source
        figure
        hold on
        
        [frequency, f_series] = spect(sensor.data(75,:),1/sensor.kgrid.dt);
        plot(frequency/1e6,f_series,'b')
        
        [frequency, f_series] = spect(sensor_data_smoothedge(75,:),1/sensor.kgrid.dt);
        plot(frequency/1e6,f_series,'g')
        
        [frequency, f_series] = spect(sensor_data_filtered(75,:),1/sensor.kgrid.dt);
        plot(frequency/1e6,f_series,'r')
        
        [frequency, f_series] = spect(simu.source.p(750,:),1/simu.kgrid.dt);
        plot(frequency/1e6,f_series,'k')
        
            title([scattering_type ' c ' num2str(c_scatt) ' rho ' num2str(rho_scatt) ',' ...
                                   ' freq ' num2str(centre_freq/1e6) ' +/- ' num2str(bandwidth/1e6)])
            xlim([0,100])
            xlabel('frequency / MHz')
            ylabel('amplitude')
            legend('sensor data - unfiltered','sensor data - smooth edge','sensor data - filtered','source')
        
        % recon with narrowband data and show image
        image_filtered = reconstruct2dUSimage(sensor_data_filtered, sensor.params, c0);
        
        figure
        imagesc(image.kgrid.x_vec*1e3,image.t_array(150:end-150)*c0*1e3,image_filtered(:,150:end-150)')
            axis image
            title([scattering_type ' c ' num2str(c_scatt) ' rho ' num2str(rho_scatt) ',' ...
                                   ' freq ' num2str(centre_freq/1e6) ' +/- ' num2str(bandwidth/1e6)])
            xlabel('x position / mm')
            ylabel('y position / mm')
            colorbar
        
        
        %% plot image
        
        figure
        imagesc(image.kgrid.x_vec*1e3,image.t_array*c0*1e3,image.data')
            axis image
            title([scattering_type ' c ' num2str(c_scatt) ' rho ' num2str(rho_scatt)])
            xlabel('x position / mm')
            ylabel('y position / mm')
            colorbar
        
        
        %% scattering distributions in hole & tmm, plot if wanted
        
        plot_toggle = false;
        
        if plot_toggle == true
            fig_distr = figure;
                title('scattering distributions')
                hold on
                xlabel('pixel intensity')
                ylabel('count')
        end
        
        [scatter_hole_mean, scatter_hole_std] = get_scattering_distr_in_hole(image, plot_toggle);
        [scatter_tmm_mean, scatter_tmm_std] = get_scattering_distr_in_tmm(image, plot_toggle);
        
        if plot_toggle == true
            legend(gca,'show')
        end
        
        
        %% image quality metrics
        
        scatSNR = scatter_tmm_mean / scatter_hole_mean;
        scatCNR = (scatter_tmm_mean - scatter_hole_mean) / (scatter_tmm_std + scatter_hole_std);
        
        disp('  scatSNR   scatCNR')
        disp([scatSNR,scatCNR])
        
        % pause
        
        
    end
end


%% compare with experimental data
% 181204 atmm with orgasol

file_dir = '../data/imagingUS/';
file_name = '181204/atmm_orgasol1_BK31[CNT]@0nm_t0[0]_dx[100�m]_dy[100�m]_dt[8ns]_03s08m21h_04-12-18_avg1_2D_raw.SGL';
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

[frequency, f_series] = spect(squeeze(exp.data(75,75,:)),1/exp.params.dt); % ,'Window','Tukey'
figure, plot(frequency/1e6,f_series,'b')
    xlim([0,100])
    ylim([0,0.05])
    xlabel('frequency / MHz')
    ylabel('amplitude')


%% LOCAL FUNCTIONS

function ROI = get_discROI_at_hole_in_image(image, fractional_radius)

    hole_x = 0;
    hole_y = 2.5e-3;
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
    
    binwidth   = 10;    
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


