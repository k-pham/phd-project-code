function obj = correcting_t0(obj, correction)
    if correction > 0
        obj.sensor_data = cat(3, zeros(obj.Nx,obj.Ny,correction), obj.sensor_data);
    elseif correction < 0
        obj.sensor_data = obj.sensor_data(:,:,-correction+1:end);
    end
    obj.Nt = obj.Nt + correction;
end