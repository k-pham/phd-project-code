function timeseries = removeDCoffset(timeseries,DC_proportion)
% removes DC offset at the start of a timeseries signal (proportion given
% in percent)

    num_avg_samples = round(length(timeseries)*DC_proportion/100);
    DCoffset = mean(timeseries(1:num_avg_samples));
    timeseries = timeseries - DCoffset;

end