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
        function obj = load_data_and_params(obj)
            [obj.sensor_data, obj.params] = loadSGL([obj.file_dir obj.file_name]);
            
            obj.Nx = obj.params.Nx;
            obj.Ny = obj.params.Ny;
            obj.Nt = obj.params.Nt;
            obj.dx = obj.params.dx;
            obj.dy = obj.params.dy;
            obj.dt = obj.params.dt;
            
            obj.Nt_delay = obj.trigger_delay / obj.dt;
        end
        
        % zero padding source
        function obj = zero_padding_source(obj, pads)
            obj.sensor_data = cat(3, zeros(obj.Nx,obj.Ny,pads), obj.sensor_data(:,:,pads+1:end) );
        end
        
        % zero padding delay
        function obj = zero_padding_delay(obj, pads)
            obj.sensor_data = cat(3, zeros(obj.Nx,obj.Ny,int32(pads)), obj.sensor_data );
            obj.Nt = obj.Nt + pads;
        end
        
        % t0 correction
        function obj = correcting_t0(obj, correction)
            if correction > 0
                obj.sensor_data = cat(3, zeros(obj.Nx,obj.Ny,correction), obj.sensor_data);
            elseif correction < 0
                obj.sensor_data = obj.sensor_data(:,:,-correction+1:end);
            end
            obj.Nt = obj.Nt + correction;
        end
        
        % zero padding sides
        function obj = zero_padding_sides(obj, pads)
            sensor_data_padded = zeros(obj.Nx+2*pads, obj.Ny+2*pads, obj.Nt);
            sensor_data_padded(pads+1:obj.Nx+pads, pads+1:obj.Ny+pads, :) = obj.sensor_data;
            assert(isequal( size(sensor_data_padded), size(obj.sensor_data)+[2*pads,2*pads,0] ))

            obj.sensor_data = sensor_data_padded;
            obj.Nx = obj.Nx + 2*pads;
            obj.Ny = obj.Ny + 2*pads;
        end
        
        % upsampling data *2
        function obj = upsampling_x2(obj)
            sensor_data_upsampled = zeros( 2*obj.Nx, 2*obj.Ny, obj.Nt );
            for i = 1:obj.Nx
                for j = 1:obj.Ny
                    sensor_data_upsampled(2*i-1,2*j-1,:) = obj.sensor_data(i,j,:);
                    sensor_data_upsampled(2*i-1,2*j  ,:) = obj.sensor_data(i,j,:);
                    sensor_data_upsampled(2*i  ,2*j-1,:) = obj.sensor_data(i,j,:);
                    sensor_data_upsampled(2*i  ,2*j  ,:) = obj.sensor_data(i,j,:);
                end
            end
            assert(isequal( size(sensor_data_upsampled), [2*obj.Nx, 2*obj.Ny, obj.Nt] ))
            
            obj.sensor_data = sensor_data_upsampled;
            obj.Nx = 2*obj.Nx;
            obj.Ny = 2*obj.Ny;
            obj.dx = obj.dx/2;
            obj.dy = obj.dy/2;
        end
        
        % apodising data
        function obj = apodising(obj)
            win = getWin([obj.Nx obj.Ny], 'Cosine');
            win = win + 0.5;
            obj.sensor_data = bsxfun(@times, win, obj.sensor_data);
        end

    end
    
end