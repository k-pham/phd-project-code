function obj = prepare_data(obj, varargin)

    % set usage defaults
    num_req_input_variables = 0;
    zero_pad_sides = 0;
    toUpsample = false;
    toApodise = false;
    freqfilter_params = {};

    % replace with user defined values if provided
    if nargin < num_req_input_variables
        error('Incorrect number of inputs.');
    elseif ~isempty(varargin)
        for input_index = 1:2:length(varargin)
            switch varargin{input_index}
                case 'ZeroPad'
                    zero_pad_sides = varargin{input_index + 1};
                    assert(isnumeric(zero_pad_sides),'Need number for number of pads.')
                case 'Upsample'
                    toUpsample = varargin{input_index + 1};
                case 'Apodise'
                    toApodise = varargin{input_index + 1};
                case 'FreqBandFilter'
                    freqfilter_params = varargin{input_index + 1};
                    assert(iscell(freqfilter_params),'Need cell array {centre_freq(Hz) bandwidth(Hz)}.')
                    if ~isempty(freqfilter_params)
                        assert(length(freqfilter_params)==2,'Need cell array of length 2.')
                    end
                otherwise
                    error('Unknown optional input.');
            end
        end
    end

    % zero pad source noise in the time series
    if obj.Nt_zero_pad_source
        obj.zero_padding_source(obj.Nt_zero_pad_source);
    end

    % add zero padding for delay
    if obj.Nt_delay ~= 0
        obj.zero_padding_delay(obj.Nt_delay);
    end

    % add/remove samples from sensor_data for t0 correction
    if obj.Nt_t0_correct
        obj.correcting_t0(obj.Nt_t0_correct);
    end

    % zero pad the sides
    if zero_pad_sides
        obj.zero_padding_sides(zero_pad_sides);
    end

    % upsample along space (x and y)
    if toUpsample
        obj.upsampling_x2;
    end

    % apodising data (to remove edge wave artefacts)
    if toApodise
        obj.apodising;
    end

    % frequency band filtering data (to use for frequency compounding)
    if ~isempty(freqfilter_params)
        obj.freq_filtering(freqfilter_params);
    end

end