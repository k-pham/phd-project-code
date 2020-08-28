function [frequency, f_series_avg] = freqSpecSGLavg(dataSGL,freq_sampling)
% average over frequency spectra of all time series using k-wave.spect and
% plot live time & frequency series and running average frequency spectrum

Nx = size(dataSGL,1);
Ny = size(dataSGL,2);
Nt = size(dataSGL,3);

f_series_avg = zeros(1,Nt);

figure(1)
set(gcf,'Position',[700 300 900 500])

for slice_x = 1:Nx
    for slice_y = 1:Ny
        t_series = squeeze(dataSGL(slice_x,slice_y,:));
        t_series = removeDCoffset(t_series,100);
        [frequency, f_series] = spect(t_series,freq_sampling,'Window','Tukey');
        f_series_avg = f_series_avg + f_series;

%         subplot(3,1,1)
%         plot(t_series)
%             title(['single time series @ x', num2str(slice_x), ' y', num2str(slice_y)])
%             xlabel('time / dt')
%             ylabel('signal amplitude / V')
%             axis([0,Nt,-0.3,0.3])
%         subplot(3,1,2)
%         semilogy(frequency/1e6, f_series)
%             title('single frequency spectrum')
%             xlabel('frequency / MHz')
%             ylabel('signal amplitude / V')
%         subplot(3,1,3)
%         semilogy(frequency/1e6, f_series_avg)
%             title('average frequency spectrum')
%             xlabel('frequency / MHz')
%             ylabel('signal amplitude / V')        
    end
end

% average
f_series_avg = f_series_avg / (Nx*Ny);

% normalise to max
f_series_avg = f_series_avg / max(f_series_avg);

% subplot(3,1,3)
semilogy(frequency/1e6, f_series_avg)
    title('average frequency spectrum')
    xlabel('frequency / MHz')
    ylabel('signal amplitude / V')
    %xlim([0 1])    
    %ylim([0 1])
    set(gca,'FontSize',13)
hold on

end
