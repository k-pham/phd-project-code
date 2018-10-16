function viewSGLpeakHisto(dataSGL,varargin)
% plot histogram of peaks along 3rd axis (/timeseries from SGL data)

    % set usage defaults
    num_req_input_variables = 1;
    binwidth = 0.01;
    color = 'standard';

    % replace with user defined values if provided
    if nargin < num_req_input_variables
        error('Incorrect number of inputs.');
    elseif ~isempty(varargin)
        for input_index = 1:2:length(varargin)
            switch varargin{input_index}
                case 'BinWidth'
                    binwidth = varargin{input_index + 1};
                case 'Color'
                    color = varargin{input_index +1};
                otherwise
                    error('Unknown optional input.');
            end
        end
    end

    Nx = size(dataSGL,1);
    Ny = size(dataSGL,2);

    peaks = max(dataSGL,[],3);
    
    binedges = 0:binwidth:max(max(peaks))*1.1;
    bincount = histcounts(peaks,binedges);
    bincentres = 0.5*(binedges(2:end)+binedges(1:end-1));
        
    figure(103)
    set(gcf,'Position',[40 500 600 400])
    hold on
%     histogram(peaks,'BinWidth',binwidth,'DisplayStyle','stairs')
    if strcmp(color,'standard')
        plot(bincentres,bincount)
    else
        plot(bincentres,bincount,'Color',color)
    end
        xlabel('peak signal')
        ylabel(['count (out of ' num2str(Nx*Ny) ')'])
        title('histogram of peak signal for each interrogation point')
    
end