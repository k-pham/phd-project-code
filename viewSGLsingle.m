function viewSGLsingle(file_dir,file_name,t_0,varargin)
% plot single time series from .txt file

    % set usage defaults
    num_req_input_variables = 3;
    toNormalise = false;
    toUseTimeAxis = true;

    % replace with user defined values if provided
    if nargin < num_req_input_variables
        error('Incorrect number of inputs.');
    elseif ~isempty(varargin)
        for input_index = 1:2:length(varargin)
            switch varargin{input_index}
                case 'Norm'
                    toNormalise = varargin{input_index + 1};
                case 'timeAxis'
                    toUseTimeAxis = varargin{input_index + 1};
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

    % extract time and signal from data
    time = dataSGLsingle(1,:);
    vAC  = dataSGLsingle(2,:);

    % modify time and signal to standardise (get real time and invert signal)
    time = (t_0 + time*1e-9) * 1e6;   % in us
    % vAC  = -vAC;

    % normalise signal if requested
    if toNormalise == true
        vAC = vAC / max(vAC);
    end

    % plot
    figure(1)
    set(gcf,'Position',[200 600 700 400])
    switch toUseTimeAxis
        case true
            plot(time,vAC)
        case false
            plot(vAC)
    end
    hold on
        switch toNormalise
            case true
                title({'normalised SPI time series',file_name},'Interpreter','none')
                ylim([-0.1 1.1])
            case false
                title({'SPI time series',file_name},'Interpreter','none')
%                 ylim([-0.1 1])    % NL prop range
    %             ylim([-1 1])    % laserGenUS range
        end
        xlabel('time / \mu s')
        ylabel('signal amplitude / V')
%         xlim([5.70 5.90])    % NL prop range
    %     xlim([6 16])    % laserGenUS range
    
end
