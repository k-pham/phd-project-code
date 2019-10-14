%% load data and choose profile

% [sensor_data, params] = loadSGL('D:\PROJECT\data\scannerCharac\191002\reflection_lineTarget_trolley.SGL');
% profile = sensor_data(6,:);

[sensor_data, params] = loadSGL('D:\PROJECT\data\scannerCharac\191014\reflection_lineTarget_benchTopScanner.SGL');
profile = sensor_data(11,:);

dprofiledx = diff(profile);

yvec = params.dy*(1:1:params.Ny);


%% plot profile and derivative

% figure; plot(yvec,profile)
figure; plot(yvec(1:end-1),dprofiledx)
hold on


%% define boundaries for segments

% bounds = round(60:200.4:20320);     % trolley scanner
bounds = round(180:200.4:20830);    % bench-top scanner
bounds_yvec = yvec(bounds);

plot(bounds_yvec,zeros(size(bounds_yvec)),'gx')

num_segments = length(bounds)-1;

peakPosition = zeros(1,num_segments);
spotSizeEsti = zeros(1,num_segments);


%% loop through segments and fit

for idx = 1:num_segments

segment = bounds(idx):bounds(idx+1)-1;
segment_data = dprofiledx(segment);
segment_yvec = yvec(segment);
segment_centre = segment_yvec(round(length(segment_yvec)/2));

if mean(segment_data) < 0
    startPoint = [-0.06 segment_centre 2e-5];
elseif mean(segment_data) >0
    startPoint = [0.06 segment_centre 2e-5];
end

f = fit(segment_yvec', segment_data', 'gauss1', 'Start', startPoint);

peakPosition(idx) = f.b1;
spotSizeEsti(idx) = f.c1*2*sqrt(log(2));

segment_fit =  f.a1*exp(-((segment_yvec-f.b1)/f.c1).^2);

plot(segment_yvec, segment_fit,'r--')

end


%%

figure
plot(peakPosition*1e3,spotSizeEsti*1e6)
hold on
    xlabel('y position / mm')
    ylabel('spot size estimate / \mum')
    xlim([0,21])
    ylim([0,70])








