function obj = zero_padding_sides(obj, pads)
    sensor_data_padded = zeros(obj.Nx+2*pads, obj.Ny+2*pads, obj.Nt);
    sensor_data_padded(pads+1:obj.Nx+pads, pads+1:obj.Ny+pads, :) = obj.sensor_data;
    assert(isequal( size(sensor_data_padded), size(obj.sensor_data)+[2*pads,2*pads,0] ))

    obj.sensor_data = sensor_data_padded;
    obj.Nx = obj.Nx + 2*pads;
    obj.Ny = obj.Ny + 2*pads;
end