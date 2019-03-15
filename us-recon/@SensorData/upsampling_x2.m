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