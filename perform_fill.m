base_dirr        = '/asap3/petra3/gpfs/p07/2019/data/11006746/';
name_ip          = 'ip_700_allpos_02';
pars_dirr        = [base_dirr,'processed/inputpar/'];

addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/core/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/additional/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/find_ub_matrix/']);
load([pars_dirr,name_ip,'.mat'],'save_dirr');

myfiles={'orthoiii2_700_merged_v2_punched_box92_sg3_bragg0p27_0p6_0p6.mat';};
         
mykernels=[5];
mysigma  =[0.3];

for idxfile=1:size(myfiles,1)
    loadnm=myfiles{idxfile};
    load([save_dirr,loadnm],'Int_total','z_total','idx_paf','hmin','hmax',...
                            'kmin','kmax','lmin','lmax','ub');
    Int_total_t=Int_total./z_total; Int_total_t(isnan(Int_total_t))=0;
    data1=Int_total_t; data1(idx_paf)=nan; clear Int_total_t
    n_steps=size(data1);
    qh=linspace(hmin,hmax,n_steps(1));
    qk=linspace(kmin,kmax,n_steps(2));
    ql=linspace(lmin,lmax,n_steps(3));
    
    % preparing for filling with python
    pyenv();
    py.importlib.import_module('astropy.convolution');
    [xx,yy,zz]=meshgrid((-mykernels(idxfile):mykernels(idxfile)));
    % usually a box of 5x5x5 (-2:2) voxels works well, but for bigger holes
    % one might want to increase that. For larger kernels than 11x11x11
    % (-5:5), it is better (for speed) to use the fft_convolve "engine" in
    % the below interpolate_replace_nans algorithm.
    sigma=mysigma(idxfile); % works best with small sigma (0:1]
    kernel=kernel_gauss3D(xx,yy,zz,0,0,0,sigma,sigma,sigma);

    % filling only the nans
    tic;
    data2=py.astropy.convolution.interpolate_replace_nans(py.numpy.array(data1),kernel);
    data2=double(data2);
    endtime=toc;
    disp(['astropy filling in ',num2str(endtime/60),' minutes.'])

    % saving
    Int_total_t=data2;
    Int_total_t(isnan(Int_total_t))=0;
    save_name=[loadnm(1:end-4),'_filled_sg',num2str(sigma),'_kernel',num2str(max(xx(:))),'.mat'];
    save([save_dirr,save_name],'Int_total_t','Int_total','z_total','idx_paf','hmin','hmax',...
                               'kmin','kmax','lmin','lmax','ub','-v7.3')
end
%% comparizon plotting with unfilled
%% k-l plot
n = 0; % q value along the integration axis
dlq=0.025; % half of the integration size

figure(201)
clf, dime(45,20)
subplot(1,2,1)
[dp,xp,yp,zp]=slice_cube(data1,qh,qk,ql,[-dlq dlq]+n,[kmin kmax],[lmin lmax],[1 0 0]);
imagesc(yp,xp,dp)
colormap(viridis), caxis([0 1]*0.7), colorbar
xlim([-10 10]),ylim([-10 10])
xlabel('{\it{l}} in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
ylabel('k in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
title([strrep(loadnm,'_','\_'),newline,...
      'h = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dlq),' r.l.u.'],...
      'interpreter','latex')
axis square, grid on
set(gca,'fontsize',17,'xtick',-100:10:100,'ytick',-100:5:100,'linewidth',1,...
    'xcolor','k','ycolor','k')

subplot(1,2,2)
[dp,xp,yp,zp]=slice_cube(data2,qh,qk,ql,[-dlq dlq]+n,[kmin kmax],[lmin lmax],[1 0 0]);
imagesc(yp,xp,dp)
colormap(viridis), caxis([0 1]*0.7), colorbar
xlim([-10 10]),ylim([-10 10])
xlabel('{\it{l}} in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
ylabel('k in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
title(['astropy filled',' sg = ',num2str(sigma),', kernel = ',num2str(max(xx(:))),newline,...
      'h = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dlq),' r.l.u.'],...
      'interpreter','latex')
axis square, grid on
set(gca,'fontsize',17,'xtick',-100:10:100,'ytick',-100:5:100,'linewidth',1,...
    'xcolor','k','ycolor','k')
% savepng([save_dirr,'plots/'],[save_name(1:end-4),'_kl'],500)

%% h-l plot
n = -3; % q value along the integration axis
dlq=0.025; % half of the integration size

figure(202)
clf, dime(45,20)
subplot(1,2,1)
[dp,xp,yp,zp]=slice_cube(data1,qh,qk,ql,[hmin hmax],[-dlq dlq]+n,[lmin lmax],[0 1 0]);
imagesc(yp,xp,dp)
colormap(viridis), caxis([0 1]*0.7), colorbar
xlim([-10 10]),ylim([-10 10])
xlabel('{\it{l}} in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
ylabel('h in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
title([strrep(loadnm,'_','\_'),newline,...
       'k = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dlq),' r.l.u.'],...
       'interpreter','latex')
axis square, grid on
set(gca,'fontsize',17,'xtick',-100:10:100,'ytick',-100:5:100,'linewidth',2,...
    'xcolor','k','ycolor','k')

subplot(1,2,2)
[dp,xp,yp,zp]=slice_cube(data2,qh,qk,ql,[hmin hmax],[-dlq dlq]+n,[lmin lmax],[0 1 0]);
imagesc(yp,xp,dp)
colormap(viridis), caxis([0 1]*0.7), colorbar
xlim([-10 10]),ylim([-10 10])
xlabel('{\it{l}} in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
ylabel('h in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
title(['astropy filled',' sg = ',num2str(sigma),', kernel = ',num2str(max(xx(:))),newline,...
       'k = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dlq),' r.l.u.'],...
       'interpreter','latex')
axis square, grid on
set(gca,'fontsize',17,'xtick',-100:10:100,'ytick',-100:5:100,'linewidth',2,...
    'xcolor','k','ycolor','k')
% savepng([save_dirr,'plots/'],[save_name(1:end-4),'_hl'],500)

%% h-k plot
n = -0.5; % q value along the integration axis
dlq=0.025; % half of the integration size

figure(203)
clf, dime(45,20)
subplot(1,2,1)
[dp,xp,yp,zp]=slice_cube(data1,qh,qk,ql,[hmin hmax],[kmin kmax],[-dlq dlq]+n,[0 0 1]);
imagesc(yp,xp,dp)
colormap(viridis), caxis([0 1]/2), colorbar
xlim([-10 10]),ylim([-10 10])
xlabel('k in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
ylabel('h in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
title([strrep(loadnm(1,:),'_','\_'),newline,...
      '{\it{l}} = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dlq),' r.l.u.'],...
      'interpreter','latex')
axis square, grid on
set(gca,'fontsize',17,'xtick',-100:5:100,'ytick',-100:5:100,'linewidth',1,...
    'xcolor','k','ycolor','k')

subplot(1,2,2)
[dp,xp,yp,zp]=slice_cube(data2,qh,qk,ql,[hmin hmax],[kmin kmax],[-dlq dlq]+n,[0 0 1]);
imagesc(yp,xp,dp)
colormap(viridis), caxis([0 1]/2), colorbar
xlim([-10 10]),ylim([-10 10])
xlabel('k in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
ylabel('h in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
title(['astropy filled',' sg = ',num2str(sigma),', kernel = ',num2str(max(xx(:))),newline,...
       '{\it{l}} = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dlq),' r.l.u.'],...
       'interpreter','latex')
axis square, grid on
set(gca,'fontsize',17,'xtick',-100:5:100,'ytick',-100:5:100,'linewidth',1,...
    'xcolor','k','ycolor','k')
% savepng([save_dirr,'plots/'],[save_name(1:end-4),'_hk'],500)

