function obj = zero_padding_source(obj, pads)
    obj.sensor_data = cat(3, zeros(obj.Nx,obj.Ny,pads), obj.sensor_data(:,:,pads+1:end) );
end