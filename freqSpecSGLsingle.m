function freqSpecSGLsingle(file_dir,file_name,freq_sampling,t_min,t_max,varargin)
% plot frequency spectrum of single time series using k-wave.spect

% set usage defaults
num_req_input_variables = 5;
toNormalise = false;

% replace with user defined values if provided
if nargin < num_req_input_variables
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Norm'
                toNormalise = varargin{input_index + 1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

% import data
fileID = fopen([file_dir file_name],'r');
formatSpec = '%f %f';
sizeData = [2 Inf];
dataSGLsingle = fscanf(fileID,formatSpec,sizeData);
fclose(fileID);

% extract time series from data
t_series = dataSGLsingle(2,t_min:t_max);

% remove DC base line (average of 10%-pre-peak signal)
avg_DCoffset = mean(t_series(1:round(length(t_series)/10)));
t_series = t_series - avg_DCoffset;

% filter t_series with Tukey window
win = getWin(length(t_series),'Tukey');
t_series = bsxfun(@times, win', t_series);

% FFT time series to get frequency spectrum using spect
[frequency, f_series] = spect(t_series,freq_sampling);

% normalise frequency spectrum if requested
if toNormalise == true
    f_series = f_series / max(f_series);
end

% plot
figure(101)
set(gcf,'Position',[200 20 700 500])
semilogy(frequency/1e6, f_series)
    %plot(frequency/1e6, 20*log(f_series))
    %plot(frequency/1e6, f_series)
hold on
    switch toNormalise
        case true
            title({'normalised freqSpec of SPI time series',file_name},'Interpreter','none')
            ylim([1e-3 1])
        case false
            title({'freqSpec of SPI time series',file_name},'Interpreter','none')
            ylim([1e-5 1e-1])
    end
    xlabel('frequency / MHz')
    ylabel('signal amplitude / V')
    %xlim([0 1])

end
