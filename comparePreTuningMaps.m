function comparePreTuningMaps(file_dir,file_dir1, file_dir2)
% Compares two pre-tuning maps by subtracting maps of (1) bias wavelength,
% (2) fitted sensitivity, (3) measured reflectivity at bias wavelength.
% Map [file_dir1] - [file_dir2]. Used for heatingEffects on pre-tuning.
    
    file_dir = {[file_dir file_dir1];
                [file_dir file_dir2]};
            
	ACgain = 1;

    % initialise arrays to save pre-tuning maps for both sets of data
    file_name_bias = 'R-map_pretuned_opt-parameters.txt';
    file_path = [file_dir{1,1} file_name_bias];
    data = textread(file_path);
    Nx = data(1,3); %mm
    Ny = data(1,6); %mm

    opts_map_full = zeros(Ny,Nx,2);
    optw_map_full = zeros(Ny,Nx,2);
    optr_map_full = zeros(Ny,Nx,2);
    
    for k = 1:2
        file_dir_k = file_dir{k,1};
        file_name_bias = 'R-map_pretuned_opt-parameters.txt';
        file_name_sens = 'Opt Rdash array_data.txt';
        file_name_sens_w = 'Opt Rdash array_wavelength.txt';
        file_name_refle = 'Opt R array_data.txt';
    

        file_path = [file_dir_k file_name_bias];
        data = textread(file_path);

%         string_map = strcat(path_map,filename_map);
%         data_full = textread(string_map);
% 
%         string_sens = strcat(path_map,file_name_sens);
%         Rdash = textread(string_sens);
%         [Rlen1,Rlen2] = size(Rdash);
% 
%         string_sens_w = strcat(path_map,file_name_sens_w);
%         Rdash_w = textread(string_sens_w);
% 
%         string_r = strcat(path_map,file_name_refle);
%         R = textread(string_r);

        
        xmin = data(1,1); %mm
        dx = data(1,2); %mm
        Nx = data(1,3); %mm
        ymin = data(1,4); %mm
        dy = data(1,5); %mm
        Ny = data(1,6); %mm
        dw = data(1,8); %mm
        Nw = data(1,9); %mm

        xidx = data(2:end,12) + 1;
        yidx = data(2:end,13) + 1;

        % falling edge
        optw = data(2:end,8);
        optr = data(2:end,9);
        opts = data(2:end,10)/ACgain;

        optw_map = zeros(Ny,Nx);
        opts_map = zeros(Ny,Nx);
        optr_map= zeros(Ny,Nx);
        for i = 1: Nx*Ny
            optw_map(yidx(i),xidx(i)) = optw(i);
            optw_list = unique(optw,'stable');
            opts_map(yidx(i),xidx(i)) = opts(i);
            optr_map(yidx(i),xidx(i)) = optr(i);
        end

        for j = 1:length(optw_list)
            optw_num(j) = sum(optw_map(:) == optw_list(j));
        end

        opts_map_full(:,:,k) = opts_map;
        optw_map_full(:,:,k) = optw_map;
        optr_map_full(:,:,k) = optr_map;

        opts_map_diff = squeeze(opts_map_full(:,:,2)) - squeeze(opts_map_full(:,:,1));
        optw_map_diff = squeeze(optw_map_full(:,:,2)) - squeeze(optw_map_full(:,:,1));
        optr_map_diff = squeeze(optr_map_full(:,:,2)) - squeeze(optr_map_full(:,:,1));

        figure(99)
        set(gcf,'Position',[80 350 900 600])
        subplot(3,3,1 + (k-1)*3)
        imagesc(optw_map)
        title('optw map')
        axis square
        colorbar

        subplot(3,3,2 + (k-1)*3)
        imagesc(opts_map)
        title('opts map')
        axis square
        colorbar

        subplot(3,3,3 + (k-1)*3)
        imagesc(optr_map)
        title('optr map')
        axis square
        colorbar

        subplot(3,3,1 + 3 + (k-1)*3)
        imagesc(optw_map_diff)
        title('optw map diff')
        axis square
        colorbar
%         caxis([0,1])
        caxis([0,0.5])

        subplot(3,3,2 + 3 + (k-1)*3)
        imagesc(opts_map_diff)
        title('opts map diff')
        axis square
        colorbar

        subplot(3,3,3 + 3 + (k-1)*3)
        imagesc(optr_map_diff)
        title('optr map diff')
        axis square
        colorbar
    end

    figure
    imagesc(optw_map_diff)
    title('optw map diff')
    axis square
    colorbar
    caxis([0,0.5])
    
    opts_map1 = opts_map_full(:,:,1);
    opts_map2 = opts_map_full(:,:,2);

    bins = -2*abs(mean(opts_map1(:))):2*abs(mean(opts_map1(:)))/100:2*abs(mean(opts_map1(:)));
    hist_opts1 = hist(opts_map1(:),bins);
    hist_opts2 = hist(opts_map2(:),bins);

    bins2 = -2*max(abs(optr_map_diff(:))):2*max(abs(optr_map_diff(:)))/100:2*max(abs(optr_map_diff(:)));
    hist_optr_diff = hist(opts_map_diff(:),bins2);

%     figure
%     plot(bins ,hist_opts1,'b',bins ,hist_opts2,'r')
%     title('opts hist')
%     axis tight
% 
%     figure
%     plot(bins2,hist_optr_diff)
%     title('optr diff hist')
    % axis([-1e-3 1e-3 0 max(hist_optr_diff(:))])
    

    %%
%     for i = 1:Rlen2
%         figure(1000)
%         imagesc(squeeze(map6(:,:,i)))
%         caxis([0 2])
%         title(i)
%         colorbar
%         pause(0.1)
%     end
    
    %%
%     figure(1001)
%     for i = 1:Ny
%         for j = 1:Nx
%             subplot(1,2,1)
%             plot(squeeze(map3(i,j,:)),squeeze(map1(i,j,:)),'r',squeeze(map4(i,j,:)),squeeze(map2(i,j,:)),'b')
%             title(['Rdash: Y: ' num2str(i) ', X: ' num2str(j)])
%     %         axis ([min_w max_w -4 4])
%             subplot(1,2,2)
%             plot(squeeze(map3(i,j,:)),squeeze(map5(i,j,:)),'r',squeeze(map4(i,j,:)),squeeze(map6(i,j,:)),'b')
%             title(['R: Y: ' num2str(i) ', X: ' num2str(j)])
%     %         axis ([min_w max_w 0 2])
%             pause(0.05)
%         end
%     end
    
    %%
%     binw = 0:0.01:0.5;
%     w_hist1 = hist(w_width(:,1),binw);
%     w_hist2 = hist(w_width(:,2),binw);
%     figure(1002)
%     plot(binw,w_hist1,'r',binw,w_hist2,'b')
    
end