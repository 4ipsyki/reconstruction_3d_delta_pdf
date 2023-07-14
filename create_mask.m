% before executing the for loop some parameters need to be set properly.
% for this purpose, load and plot a few selected positions and apply the mask.
% see the result with an overlay plot.
% see the commented parts for this purposes.

base_dirr        = '/asap3/petra3/gpfs/p07/2019/data/11006746/';
name_ip          = 'ip_700_allpos_01';
pars_dirr        = [base_dirr,'processed/inputpar/'];

addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/core/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/additional/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/find_ub_matrix/']);
load([pars_dirr,name_ip,'.mat']);

% mynr=1;
% img=zeros(x_size,y_size);
% for ii=1:50
%     img=img+double(get_pe_new5(dirr{mynr}, file{mynr}, ii));
% end
%%

nr_scans=size(file,1);
mask=cell(nr_scans,1);
for mynr=1:nr_scans
    % loading basic mask
    mask_temp=imread([base_dirr,'processed/reconstructions/ybco_oleh/mask_pil2M_2023/mask_pil2M.tif']);
    idxneg=img==-2;
    mask_temp(idxneg)=1;

    x=1:x_size; y=1:y_size;
    [xx, yy]=meshgrid(y,x);

    % circle at the center
    radius = 40; % radius of the beam stop in px
    mask_temp(((xx-orgy(mynr)).^2+(yy-orgx(mynr)).^2)<=radius^2)=1;

    % line for the beamstop holder stick
    mytt=30; % tickness for the line in px
    mym=(1190-orgx(5))/(1339-orgy(5));   % slope of the line
    myq=orgx(5)+orgx(mynr)-orgx(5)-mym*(orgy(5)+orgy(mynr)-orgy(5));   % f(0) for the line
    % mym and myq are base on the line coordinates found in fifth run
    for myn=y
        for myj=x
            if (myj-mym*myn<=myq+mytt) && (myj-mym*myn>=myq-mytt) && (myn >= orgy(mynr))
                mask_temp(myj,myn)=1;
            end
        end
    end

    % chamber-border mask
    radius = 2060;
    % in this case the center of the chamber was not centerred with the sample
    % the reference run is number one
    cx=1600+orgy(mynr)-orgy(1);
    cy=x_size-1908+orgx(mynr)-orgx(1);
    mask_temp(((xx-cx).^2+(yy-cy).^2)>=radius^2)=1;

    % % border of the image
    % mybrd=1;
    % mask(x<=mybrd,:)=1;
    % mask(:,y<=mybrd)=1;
    % mask(x>=x_size-mybrd,:)=1;
    % mask(:,y>=y_size-mybrd)=1;
    % %

    % plot
    % figure(2)
    % clf
    % imagesc(img),shading flat
    % colormap jet
    % caxis([0 10])
    % set(gca,'layer','bottom')
    % hold on
    % pltmask=mask_temp; pltmask(~mask_temp)=nan;
    % pcolor(pltmask), shading flat
    % hold off
    %%
    mask{mynr}=boolean(mask_temp);
end
save([save_dirr,'masks_700.mat'],'mask')

