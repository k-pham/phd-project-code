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