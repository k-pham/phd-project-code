function [frequency, f_series_avg] = freqSpecSGLavg2(dataSGL,freq_sampling)
% frequency spectrum of average over all time series using k-wave.spect
% DOESN'T WORK - GETS RID OF ALL SIGNAL USE freqSpecSGLavg INSTEAD

% average time series in space
t_series_avg = squeeze(mean(mean(dataSGL,1),2));

figure
plot(t_series_avg)

% remove DC signal
% t_series_avg = removeDCoffset(t_series_avg,100);

% take fourier transform
[frequency, f_series_avg] = spect(t_series_avg,freq_sampling,'Window','Tukey');

% normalise to max
f_series_avg = f_series_avg / max(f_series_avg);

figure(7)
semilogy(frequency/1e6, f_series_avg)
    title('average frequency spectrum')
    xlabel('frequency / MHz')
    ylabel('signal amplitude / V')
    set(gca,'FontSize',13)
hold on

end
