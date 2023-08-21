base_dirr        = '/asap3/petra3/gpfs/p07/2019/data/11006746/';
name_ip          = 'ip_700_allpos_02';
pars_dirr        = [base_dirr,'processed/inputpar/'];

addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/core/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/additional/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/find_ub_matrix/']);
load([pars_dirr,name_ip,'.mat']);

clear mask; load([save_dirr,'masks_700.mat']);
% clear mask; mask=imread([base_dirr,'mask_pil2M.tif']);

[xx,yy]=meshgrid(1:y_size,1:x_size);

% parameters for the automatic the automatic masking of reagions containing peaks
radius0=10;   % marking-circle radius for high Q and very small Q
radius1=100;  % masking-circle radius for Q = d0
d0=200;       % distance from center (Q) in px for the first strong peak
d1=700;       % distance from center in px for the high Q region

nr_scans=size(file,1);
abs_corr=cell(nr_scans,1);
for nr=1:nr_scans
    nr_frames=size(frame_nr{nr},2);
    % suming all framnes to find "empty" regions
    img=zeros(x_size,y_size);
    parfor ii=1:nr_frames
        if data_type==1
            img = img+single(get_pe_new5(dirr{nr}, file{nr}, ii));
        elseif data_type==2
            img = img+permute(double(read_cbf([dirr{nr}, file{nr}, sprintf('%05d.cbf',ii)]).data),[2 1]);
        else
            disp("Data type not specified properly.")
            break
        end
    end
    imgtot=img;
    mask0=mask{nr};
    
    figure(51)
    clf,imagesc(img), caxis([0 1e4])
    % isolating the "empty" regions
    peaks=FastPeakFind(img, 500);
    hold on
    plot(peaks(1:2:end),peaks(2:2:end),'r+')
    hold off
    
    myslope=(radius1-radius0)/(d0-d1);
    myq=radius1-myslope*d0;
    img1=img;
    mask1=mask0;
    for ii=1:2:length(peaks)
        mydist=sqrt((orgy(nr)-peaks(ii))^2+(orgx(nr)-peaks(ii+1))^2);
        if mydist < d0 || mydist > d1
            curr_radius = radius0;
        else
            curr_radius = myslope*mydist + myq;
        end
        mask1(((xx-peaks(ii)).^2+(yy-peaks(ii+1)).^2)<=curr_radius^2)=1;
    end
    img1(mask1)=0;
    figure(52)
    % clf,imagesc(mask1), caxis([0 1])
    clf,imagesc(img1), caxis([0 2000])
    
    % extracting the background intensity variation
    int_roi=zeros(nr_frames,1);
    parfor ii=1:nr_frames
        img=single(get_pe_new5(dirr{nr}, file{nr}, ii));
        disp(num2str(ii))
        bck=img(~mask1);
        int_roi(ii,1)=mean(bck(:));
    end

    int_roi_tmp=int_roi;
    for box_size=[1000,500,100,10] % average-moving box sizes
        box_frac=1; % keep at 1! - fraction on the box to be used as the step
        box_step=round(box_size*box_frac);
        box_perc=0.1; % threshold in fraction of median to be used as outlier removal

        box_tmp=nr_frames/box_step;
        if abs(round(box_tmp)-box_tmp)>0
            box_iter=floor(box_tmp)+1;
        else
            box_iter=box_tmp;
        end


        for ii=1:box_iter
            disp([1,box_size]+box_step*(ii-1))
            if (box_size+box_step*(ii-1))>nr_frames
                box_med_tmp=int_roi_tmp(box_step*(ii-1):end);
            else
                box_med_tmp=int_roi_tmp((1:box_size)+box_step*(ii-1));
            end
            idx_nan=isnan(box_med_tmp);
            box_idx_tmp=find((box_med_tmp(~idx_nan)/median(box_med_tmp(~idx_nan))-1)>box_perc);
            int_roi_tmp(box_idx_tmp+box_step*(ii-1))=nan;
        end
        tmp=int_roi_tmp(1:200);
        int_roi_tmp(1)=median(tmp(~isnan(tmp)));
        idx_nan=isnan(int_roi_tmp);
        int_roi_tmp=pchip(find(~idx_nan),int_roi_tmp(~idx_nan),(1:nr_frames));
    end
    tmp=int_roi_tmp(1:200);
    int_roi_tmp(1)=median(tmp(~isnan(tmp)));
    idx_nan=isnan(int_roi_tmp);
    int_spline=pchip(find(~idx_nan),int_roi_tmp(~idx_nan),(1:nr_frames));

    figure(53)
    clf, dime(15,15)
    hold on
    plot(int_roi,'-sb')
    plot(int_roi_tmp,'-r','linewidth',1.5)
    plot(smooth(int_spline),'-g','linewidth',1.5)
    hold off
    xlim([1 nr_frames]),ylim([0.0 0.5])

    % saving the result
    abs_corr{nr}=int_spline/0.05;
end
%% renormalization
% the scaling factor for the structure-factor variation is already included
% in the raw intensity. So we need only to scale each frame to the sample
% volume fraction and beam-flux variation. We consider that the conditions
% (sample volume and flux) are similar for the first 100 frames.
for ii=1:nr_scans
    scl_fkt=mean(abs_corr{ii}(1:100));
    disp(scl_fkt)
    abs_corr{ii}=abs_corr{ii}/scl_fkt;
end
%% saving result
save([save_dirr,'abs_corr_700.mat'],'abs_corr')

%%
figure(54)
clf, hold on
for ii=1:nr_scans
    plot(abs_corr{ii})
end
box on,grid on,xlim([1 nr_frames])
title('background variation for all positions')

