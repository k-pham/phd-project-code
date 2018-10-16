function freqSpecSGL(dataSGL,freq_sampling,slice_x,slice_y,t_min,t_max)
% plot frequency spectrum of selected time series using k-wave.spect

t_series = squeeze(dataSGL(slice_x,slice_y,t_min:t_max));
[frequency, f_series] = spect(t_series,freq_sampling); %,'Window','Cosine');

figure
%plot(frequency/1e6, 20*log(f_series))
%plot(frequency/1e6, f_series)
semilogy(frequency/1e6, f_series)
    title('frequency spectrum of selected time series')
    xlabel('frequency / MHz')
    ylabel('signal amplitude / V')
    %xlim([0 1])    
    %ylim([0 1])

end
