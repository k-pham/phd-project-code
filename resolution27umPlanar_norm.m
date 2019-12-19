%% find c0 and t0 given some more or less arbitrary norm

file_dir = 'D:\PROJECT\figures\_Matlab figs\USimaging\191126 resolution27umPlanar BK31[CNT] trolley scrambled fibre centralised parallel phantom\';
file_dir = [file_dir 'x-05 c0_var t0_var\'];

c0_sample = 1484:1:1494; % in m/s
t0_sample = -18:1:-12;

dist_lims = 22;%2:0.1:10; % in mm

c0_minnorms = zeros(length(t0_sample),length(dist_lims));
% t0_minnorms = zeros(1,length(dist_lims));

for idx_t = 1:length(t0_sample)

    for idx_dist = 1:length(dist_lims)

        dist_lim = dist_lims(idx_dist);

        Norm = zeros(length(t0_sample),length(c0_sample));
        Norm_std = zeros(length(t0_sample),length(c0_sample));

        for idx_c = 1:length(c0_sample)

            c0 = c0_sample(idx_c);
            t0 = t0_sample(idx_t);

            file_name = ['peaksInfoAll_x-05_c' sprintf('%0.1f',c0) '_t0_' num2str(t0) '.mat'];
            data = load([file_dir file_name]);
            peaksInfoAll = data.peaksInfoAll;

            peaksPosX    = peaksInfoAll(2,:)*1e3; % in mm
            peaksPosZ    = peaksInfoAll(3,:)*1e3; % in mm
            peaksResoLat = peaksInfoAll(4,:)*1e6; % in um

            includeInNorm = rectangularNorm2(peaksPosX, peaksPosZ, dist_lim);

            peaksResoLatInNorm = peaksResoLat(includeInNorm);

            % assert(length(peaksResoLatInNorm)==104);
            Norm(idx_t,idx_c) = mean(peaksResoLatInNorm);
            Norm_std(idx_t,idx_c) = std(peaksResoLatInNorm);

        end

        figure
        plot(c0_sample,Norm(idx_t,:))
        hold on
    %     plot(c0_sample, Norm - Norm_std, '--')
    %     plot(c0_sample, Norm + Norm_std, '--')
        xlabel('c0 [m/s]')
        ylabel('t0 correction [dt]')
    %     ylabel('avg. lat. resolution within norm-radius [\mum]')
        title(['norm radius = ' num2str(dist_lim) ' mm'])
        set(gca,'FontSize',13)
    %     axis([1460,1500,40,150])
        drawnow
        pause(0.1)

%         minnorm = min(Norm,[],'all');
%         [idx_minnorm_t, idx_minnorm_c] = find(Norm==minnorm);
        [minnorm, idx_minnorm_c] = min(Norm(idx_t,:));
        c0_minnorms(idx_t,idx_dist) = c0_sample(idx_minnorm_c);
    %     t0_minnorms(idx_dist) = t0_sample(idx_minnorm_t);

    end
    
figure
% yyaxis left
plot(dist_lims, c0_minnorms(idx_t,:), 'x')
xlabel('norm length scale [mm]')
ylabel('c0 with minimum norm [m/s]')
% yyaxis right
% plot(dist_lims, t0_minnorms, 'o')
% ylabel('t0 correction with minimum norm [dt]')
title('best c0 for different norm radii')
set(gca,'FontSize',13)

end

figure
imagesc(dist_lims, t0_sample, c0_minnorms)

function includeInNorm = circularNorm(peaksPosX, peaksPosZ, radius)

    distance = sqrt(peaksPosX.^2+peaksPosZ.^2); % in mm

    includeInNorm = distance < radius; % in mm

end

function includeInNorm = squareNorm(peaksPosX, peaksPosZ, sidelength)

    withinX = abs(peaksPosX) < sidelength/2;
    withinZ = peaksPosZ < sidelength;
    
    includeInNorm = withinX & withinZ;

end

function includeInNorm = rectangularNorm(peaksPosX, peaksPosZ, sidelength)

    withinX = abs(peaksPosX) < sidelength;
    withinZ = peaksPosZ < sidelength;
    
    includeInNorm = withinX & withinZ;

end

function includeInNorm = rectangularNorm2(peaksPosX, peaksPosZ, sidelength)

    withinX = abs(peaksPosX) < sidelength/2;
    withinZ = peaksPosZ < sidelength*2;
    
    includeInNorm = withinX & withinZ;

end

