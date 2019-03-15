function obj = zero_padding_delay(obj, pads)
    obj.sensor_data = cat(3, zeros(obj.Nx,obj.Ny,int32(pads)), obj.sensor_data );
    obj.Nt = obj.Nt + pads;
end