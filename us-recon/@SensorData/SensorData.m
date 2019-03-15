classdef SensorData < handle
    
    properties
        file_dir
        file_name
        sensor_data
        params
        Nx
        Ny
        Nt
        dx
        dy
        dt
        trigger_delay
        Nt_delay
        Nt_zero_pad_source
        Nt_t0_correct
    end
    
    methods
        
        % constructor
        function obj = SensorData(file_dir, file_name, trigger_delay, Nt_zero_pad_source, Nt_t0_correct)
            obj.file_dir  = file_dir;
            obj.file_name = file_name;
            
            obj.trigger_delay = trigger_delay;
            obj.Nt_zero_pad_source = Nt_zero_pad_source;
            obj.Nt_t0_correct = Nt_t0_correct;
        end
        
        % load data and params
        obj = load_data_and_params(obj)
        
        % zero padding source
        obj = zero_padding_source(obj, pads)
        
        % zero padding delay
        obj = zero_padding_delay(obj, pads)
        
        % t0 correction
        obj = correcting_t0(obj, correction)
        
        % zero padding sides
        obj = zero_padding_sides(obj, pads)
        
        % upsampling data *2
        obj = upsampling_x2(obj)
        
        % apodising data
        obj = apodising(obj)

        % frequency bandpass filtering data
        obj = freq_filtering(obj, freqfilter_params)
        
        % prepare data maybe using all processing methods
        obj = prepare_data(obj, varargin)
        
    end
    
end