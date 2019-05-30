function freqSpecSGLsingle(file_dir,file_name,freq_sampling,t_min,t_max,varargin)
% plot frequency spectrum of single time series using k-wave.spect

num_avg = 12;
min_length = 2001;

% set usage defaults
num_req_input_variables = 5;
toRemoveDC = true;
toApplyTukey = true;
toCorrect4PD = true;
toNormaliseT = false;
toNormaliseF = false;
linecolour = 'b';

% replace with user defined values if provided
if nargin < num_req_input_variables
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'removeDC'
                toRemoveDC = varargin{input_index + 1};
            case 'applyTukey'
                toApplyTukey = varargin{input_index + 1};
            case 'correct4PD'
                toCorrect4PD = varargin{input_index + 1};
            case 'NormTime'
                toNormaliseT = varargin{input_index + 1};
            case 'NormFreq'
                toNormaliseF = varargin{input_index + 1};
            case 'LineColour'
                linecolour = varargin{input_index + 1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

%% import data

file_type = file_name(end-2:end);

switch file_type
    
    case 'txt'
        fileID = fopen([file_dir file_name],'r');
        formatSpec = '%f %f';
        sizeData = [2 Inf];
        dataSGLsingle = fscanf(fileID,formatSpec,sizeData);
        fclose(fileID);

        t_series = dataSGLsingle(2,t_min:t_max);
        time = dataSGLsingle(1,:) * 1e-3; % from ns to us
        
    case 'SGL'
        [dataSGL, params] = loadSGL([file_dir file_name]);
        
        dataSGL_reshaped = reshape(dataSGL, [ size(dataSGL,1)*size(dataSGL,2) size(dataSGL,3) ]);
        dataSGLsingle = mean(dataSGL_reshaped(1:num_avg, :), 1);
        assert(length(dataSGLsingle) == size(dataSGL,3));
        
        t_series = dataSGLsingle(t_min:t_max);
        time = (params.t0 + (1:params.Nt) * params.dt) * 1e6; % from s to us
        
end


%% remove DC base line (average of 10%-pre-peak signal)

if toRemoveDC == true
    t_series = removeDCoffset(t_series,10);
end


%% filter t_series with Tukey window

if toApplyTukey == true
    tukey_win = getWin(length(t_series),'Tukey');
    t_series = bsxfun(@times, tukey_win', t_series);
end


%% padding of t_series and tukey window (to ensure same size freq bins)

if length(t_series) < min_length
    num_pad_samples_front = round((min_length - length(t_series))/2);
    num_pad_samples_back = min_length - length(t_series) - num_pad_samples_front;
    
    t_series  = [ zeros(1,num_pad_samples_front) t_series   zeros(1,num_pad_samples_back) ];
	tukey_win = [ zeros(1,num_pad_samples_front) tukey_win' zeros(1,num_pad_samples_back) ];
    time = time( t_min-num_pad_samples_front : t_max+num_pad_samples_back ) * 1e-3;
end


%% FFT time series to get frequency spectrum using spect

[frequency, f_series] = spect(t_series,freq_sampling);
f_tukey = spect(tukey_win,freq_sampling);


%% correct for PD response

if toCorrect4PD == true
    % import transimpedance (V/A) gain from Edward's excel sheet
    PDresponse = xlsread('../data/Log of Photodetectors.xlsx','PD230MHz-45kOhm#01','A1:B2402');
    
    PD_freq = PDresponse(:,1) * 1e6; % convert MHz to Hz
    PD_gain = PDresponse(:,2);
	
    PD_gain_resample = interp1(PD_freq,PD_gain,frequency);
    PD_gain_resample(isnan(PD_gain_resample)) = min(PD_gain_resample);
    
    f_series = f_series ./ PD_gain_resample;
end


%% normalise time series and/or frequency spectrum if requested
if toNormaliseT == true
    t_series = t_series / max(t_series);
end

if toNormaliseF == true
    f_series = f_series / max(f_series);
end


%% plot time series

figure(101)
set(gcf,'Position',[100 500 700 400])
plot(time,t_series,linecolour)
hold on
    xlabel('time / /mu s')
    ylabel('signal amplitude / V')


%% plot frequency spectrum

figure(102)
set(gcf,'Position',[100 20 700 400])
semilogy(frequency/1e6, f_series,linecolour)
    %plot(frequency/1e6, 20*log(f_series))
    %plot(frequency/1e6, f_series)
hold on
    switch toNormaliseF
        case true
            title({'normalised freqSpec of SPI time series',file_name},'Interpreter','none')
            ylim([1e-4 1])
            xlim([0,400])
        case false
            title({'freqSpec of SPI time series',file_name},'Interpreter','none')
%             ylim([1e-5 1e-1])
    end
    xlabel('frequency / MHz')
    ylabel('signal amplitude / V')


%% plot with part of time series used, Tukey filter, freq spectra of series & filter
figure(103)
set(gcf,'Position',[200 300 1000 600])
subplot(2,2,1)
hold on
plot(time,t_series,linecolour)
    title('windowed & filtered time series')
    xlabel('time [us]')
    ylabel('signal [mV]')
%     axis([3,3.5,-100,1000])
subplot(2,2,2)
hold on
plot(time,tukey_win,linecolour)
    title('Tukey window filter')
    xlabel('time [us]')
    ylabel('signal [mV]')
	ylim([-0.1,1.1])
subplot(2,2,3)
semilogy(frequency/1e6, f_series,linecolour)
hold on
    title('frequency spectrum')
    xlabel('frequency [MHz]')
    ylabel('signal amplitude')
    xlim([0,400])
%     ylim([1e-4 1])
subplot(2,2,4)
semilogy(frequency/1e6, f_tukey,linecolour)
hold on
    title('freq spec of tukey')
    xlabel('frequency [MHz]')
    ylabel('signal amplitude')
    xlim([0,100])
    ylim([1e-5,1])

end
