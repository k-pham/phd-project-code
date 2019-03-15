function obj = apodising(obj)
    win = getWin([obj.Nx obj.Ny], 'Cosine');
    win = win + 0.5;
    obj.sensor_data = bsxfun(@times, win, obj.sensor_data);
end