function comparePeakSGL(dataSGL_1, dataSGL_2)
% pick out peak (positive maximum) of all time series and compare for 1/cold
% and 2/hot plot in scatter plot showing peak at each interrogation point

    % separate peaks
    peaks_1 = max(dataSGL_1,[],3);
    peaks_2 = max(dataSGL_2,[],3);
    
    Nx = size(dataSGL_1,1);
    Ny = size(dataSGL_1,2);
    
    % compute sensitivity loss as ratio of peaks
    sensitivity_loss = peaks_2 ./ peaks_1;
    sens_loss = reshape(sensitivity_loss,[1,Nx*Ny]);
    sens_loss_mean = mean(sens_loss);
    sens_loss_std  = std(sens_loss);
    
    % surf plot of sensitivity loss over interrogation points
%     figure
%     surf(sensitivity_loss)
%         zlim([0,2])
%         caxis([0,1])
    
    % histograms
    figure
    set(gcf,'Position',[100 400 600 600])
    % histomgram of sensitivity/peak pressures for both cold/hot overlaid
    subplot(2,1,1)
    hold on
    binwidth = 0.01;
    histogram(peaks_1,'BinWidth',binwidth)
    histogram(peaks_2,'BinWidth',binwidth)
        xlabel('peak signal')
        ylabel(['count (out of ' num2str(Nx*Ny) ')'])
        title('histogram of peak signal for each interrogation point')
        legend('peaks_{data1}','peaks_{data2}')
    hold off
    % histogram of sensitivity loss over all interrogation points
    subplot(2,1,2)
    histogram(sensitivity_loss,'BinWidth',0.02,'FaceColor',1/255*[0,0,0])
        hold on
        plot([1,1],[0,2000],'r-')
        axis([0,2,0,1500])
        xlabel('sensitivity loss = peaks_{data2} / peaks_{data1}')
        ylabel(['count (out of ' num2str(Nx*Ny) ')'])
        title('histogram of sensitivity loss for each interrogation point')
    
    % scatter plot of peaks for each interrogation point
    figure
    set(gcf,'Position',[300 0 600 500])
    peak_1_vec = reshape(peaks_1,[1,Nx*Ny]);
    peak_2_vec = reshape(peaks_2,[1,Nx*Ny]);
    scatter(peak_1_vec,peak_2_vec,'.','MarkerEdgeColor',1/255*[38,38,38])
        xlabel('peaks_{data1}')
        ylabel('peaks_{data2}')
        title({ 'correlation between peak_{1} and peak_{2} for each interrogation point', ...
%         title({ 'correlation between peak_{hot 20Hz 2} and peak_{hot 20kHz 3} for each interrogation point', ...
                ['sensitivity loss = ' num2str(sens_loss_mean) '+/-' num2str(sens_loss_std)]    })
%         xlabel('peaks_{heat20Hz}')
%         ylabel('peaks_{heat20Hz no UST}')
        axis equal
        axis([0,0.2,0,0.2])
        hold on
        plot([0,0.2],[0,0.2],'r-')
        
end