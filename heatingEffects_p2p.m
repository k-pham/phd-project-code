function heatingEffects_p2p(file_dir,file_name,dt_samples,num_samples,frac_presamples,varargin)

    % set usage defaults
    num_req_input_variables = 5;
    toSubtractBaseline = false;

    % replace with user defined values if provided
    if nargin < num_req_input_variables
        error('Incorrect number of inputs.');
    elseif ~isempty(varargin)
        for input_index = 1:2:length(varargin)
            switch varargin{input_index}
                case 'subtractBaseline'
                    toSubtractBaseline = varargin{input_index + 1};
                otherwise
                    error('Unknown optional input.');
            end
        end
    end

    fileID = fopen([file_dir file_name],'r');
    formatSpec = '%f %f %f';
    sizeData = [3 Inf];
    data_p2p = fscanf(fileID,formatSpec,sizeData);
    fclose(fileID);
    data_p2p = data_p2p';

    time = data_p2p(:,1);
    vAC  = data_p2p(:,2);
    vDC  = data_p2p(:,3);

    % make real time array in ms and setting t_trigger = 0
    realtime = linspace(dt_samples,dt_samples*num_samples,num_samples); %in seconds
    realtime = realtime*1e3;    % in ms
    realtime = realtime - realtime(end)*frac_presamples;    % subtract presamples

    % subtract baseline DC (first 10% of signal/before t = 0) from DC
    if toSubtractBaseline == true
        vDC_base = mean(vDC(realtime < 0));
        vDC = vDC - vDC_base;
    end
    
    figure(100)
    hold on
    plot(realtime,vDC(1:num_samples)) % ,'r', 'Color',[0 0.4470 0.7410]
        title(file_name,'Interpreter','none')
        xlabel('time [ms]')
        switch toSubtractBaseline
            case true
                ylabel('vDC - vDC_{baseline}')
            otherwise
                ylabel('vDC')
        end
        xlim([realtime(1),realtime(end)])
        ylim([-0.5,1])

end