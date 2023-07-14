base_dirr        = '/asap3/petra3/gpfs/p07/2019/data/11006746/';
name_ip          = 'ip_700_allpos_02';
pars_dirr        = [base_dirr,'processed/inputpar/'];

addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/core/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/additional/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/find_ub_matrix/']);
load([pars_dirr,name_ip,'.mat']);

marcopos=round(nr/4)+1;
% loading the reconstruction
loadnm=['orthoiii2_700_merged_v2.mat'];
save_name = [save_dirr,loadnm];
load(save_name,'Int_total_t','hmin','hmax','kmin','kmax','lmin','lmax');

Int_total_t(isnan(Int_total_t))=0;
n_steps=size(Int_total_t);
qh=linspace(hmin,hmax,n_steps(1));
qk=linspace(kmin,kmax,n_steps(2));
ql=linspace(lmin,lmax,n_steps(3));

%% k-l plot
n = 0; % q value along the integration axis
dlq=0.05; % half of the integration size

figure(101)
clf, dime(20,20)
[dp,xp,yp,zp]=slice_cube(Int_total_t,qh,qk,ql,[-dlq dlq]+n,[kmin kmax],[lmin lmax],[1 0 0]);
imagesc(yp,xp,dp)
colormap(viridis), caxis([0 0.8]), colorbar
xlim([lmin lmax]), ylim([kmin kmax])
xlabel('{\it{l}} in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
ylabel('k in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
title([strrep(loadnm,'_','\_'),newline,...
      'h = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dlq),' r.l.u.'],...
      'interpreter','latex')
axis square, grid on
set(gca,'fontsize',17,'xtick',-100:2:100,'ytick',-100:2:100,'linewidth',1,...
    'xcolor','k','ycolor','k')
%savepng([save_dirr,'plots/'],[loadnm,'_kl'],500)

%% h-l plot
n = 0; % q value along the integration axis
dlq=0.05; % half of the integration size

figure(102)
clf, dime(20,20)
[dp,xp,yp,zp]=slice_cube(Int_total_t,qh,qk,ql,[hmin hmax],[-dlq dlq]+n,[lmin lmax],[0 1 0]);
imagesc(yp,xp,dp)
colormap(viridis), caxis([0 0.8]), colorbar
xlim([lmin lmax]), ylim([hmin hmax])
xlabel('{\it{l}} in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
ylabel('h in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
title([strrep(loadnm,'_','\_'),newline,...
       'k = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dlq),' r.l.u.'],...
       'interpreter','latex')
axis square, grid on
set(gca,'fontsize',17,'xtick',-100:2:100,'ytick',-100:2:100,'linewidth',2,...
    'xcolor','k','ycolor','k')
%savepng([save_dirr,'plots/'],[loadnm,'_hl'],500)

%% h-k plot
n = 0; % q value along the integration axis
dlq=0.05; % half of the integration size

figure(103)
clf, dime(20,20)
[dp,xp,yp,zp]=slice_cube(Int_total_t,qh,qk,ql,[hmin hmax],[kmin kmax],[-dlq dlq]+n,[0 0 1]);
imagesc(yp,xp,dp)
colormap(viridis), caxis([0 0.8]), colorbar
xlim([kmin kmax]), ylim([hmin hmax])
xlabel('k in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
ylabel('h in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
title([strrep(loadnm,'_','\_'),newline,...
      '{\it{l}} = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dlq),' r.l.u.'],...
      'interpreter','latex')
axis square, grid on
set(gca,'fontsize',17,'xtick',-100:2:100,'ytick',-100:2:100,'linewidth',1,...
    'xcolor','k','ycolor','k')
%savepng([save_dirr,'plots/'],[loadnm,'_hk'],500)

