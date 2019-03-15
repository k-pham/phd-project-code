function obj = freq_filtering(obj, freqfilter_params)
    disp('Frequency bandpass filtering ...'),
    tic

    [centre_freq, bandwidth] = freqfilter_params{:};
    bandwidth_pc = bandwidth / centre_freq * 100;       % convert bandwidth to % of centre frequency for gaussianFilter
    %sensor_data_filtered = zeros(size(obj.sensor_data));

    for i = 1:obj.Nx
        obj.sensor_data(i,:,:) = gaussianFilter(squeeze(obj.sensor_data(i,:,:)),1/obj.dt,centre_freq,bandwidth_pc);
    end

    disp(['  completed in ' scaleTime(toc)]);
end