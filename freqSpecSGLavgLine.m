function [frequency, f_series_avg] = freqSpecSGLavgLine(dataSGL,params)
% average over frequency spectra of all time series using k-wave.spect and
% plot live time & frequency series and running average frequency spectrum

Nx = params.Nx;     % size(dataSGL,1);
Ny = params.Ny;     % size(dataSGL,1);
Nt = size(dataSGL,2); % params.Nt;
freq_sampling = 1/params.dt;

if min(Nx,Ny) ~= 1
    error('Data not from a line scan.');
end

% f_series_avg = zeros(1,Nt/2+1);

% ---
f_series_avg = zeros(1,21);
% ---

for slice = 1:max(Nx,Ny)
    
    t_series = dataSGL(slice,:);
        
    % ---
    [~ , maxIndex] = max(t_series);
    t_series = t_series(:,maxIndex-20:maxIndex+20);
    % ---
    
    [frequency, f_series] = spect(t_series,freq_sampling); %,'Window','Cosine');
    f_series_avg = f_series_avg + f_series;

    figure(7)
    set(gcf,'Position',[70 180 700 800])
    subplot(3,1,1)
    plot(t_series)
        title(['single time series @', num2str(slice)])
        xlabel('time / dt')
        ylabel('signal amplitude / V')
%         axis([0,Nt,-0.3,0.3])
    subplot(3,1,2)
    semilogy(frequency/1e6, f_series)
        title('single frequency spectrum')
        xlabel('frequency / MHz')
        ylabel('signal amplitude / V')
    subplot(3,1,3)
    semilogy(frequency/1e6, f_series_avg)
        title('average frequency spectrum')
        xlabel('frequency / MHz')
        ylabel('signal amplitude / V')        
end

f_series_avg = f_series_avg / max(Nx,Ny);

subplot(3,1,3)
semilogy(frequency/1e6, f_series_avg)
    title('average frequency spectrum')
    xlabel('frequency / MHz')
    ylabel('signal amplitude / V')
    %xlim([0 1])    
    %ylim([0 1])

end
