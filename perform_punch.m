base_dirr        = '/asap3/petra3/gpfs/p07/2019/data/11006746/';
name_ip          = 'ip_700_allpos_02';
pars_dirr        = [base_dirr,'processed/inputpar/'];

addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/core/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/additional/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/find_ub_matrix/']);
load([pars_dirr,name_ip,'.mat'],'save_dirr');

% loading the reconstruction
loadnm='orthoiii2_700_merged_v2.mat';
save_name = [save_dirr,loadnm];
load(save_name,'Int_total_t','hmin','hmax','kmin','kmax','lmin','lmax','ub');

Int_total_t(isnan(Int_total_t))=0;
data1=Int_total_t; clear Int_total_t
n_steps=size(data1);
qh=linspace(hmin,hmax,n_steps(1));
qk=linspace(kmin,kmax,n_steps(2));
ql=linspace(lmin,lmax,n_steps(3));

%% punching
data2=data1;

box_size = 92;              % moving-box size NxNxN
thr_sigm =  3;              % threshold in multiples of sigma of the data inside the moving box
fill_flg =  -1;             % statistical-flat-fill option (see paf_karen.m for more explanation)
delta_qB = [0.27 0.6 0.6];  % radius of the region around Bragg peaks (in r.l.u.) to be punched
delta_sh =  1;              % shape of the region around the Bragg: 1 - elispsoid; 2 - cube
[data2,idx_paf]=paf_karen(data2,box_size,thr_sigm,fill_flg,[qh;qk;ql],delta_qB,delta_sh);

%% saving the punched data
Int_total_t=data2; Int_total_t(isnan(Int_total_t))=0;
file_name=[loadnm(1:end-4),'_punched_box',num2str(box_size),...
                           '_sg',num2str(thr_sigm),...
                           '_bragg',strrep(strrep(num2str(delta_qB),'.','p'),'         ','_'),'         ','_'),...
                           '.mat'];
save([save_dirr,file_name],'Int_total_t','idx_paf','hmin','hmax',...
                           'kmin','kmax','lmin','lmax','ub','-v7.3');

%% comparison plotting with unpunched
%% k-l plot
n = -0; % q value along the integration axis
dlq=0.025; % half of the integration size

figure(201)
clf, dime(45,20)
subplot(1,2,1)
[dp,xp,yp,zp]=slice_cube(data1,qh,qk,ql,[-dlq dlq]+n,[kmin kmax],[lmin lmax],[1 0 0]);
imagesc(yp,xp,dp)
colormap(viridis), caxis([0 1]/2), colorbar
xlim([-20 20]./2),ylim([-10 10]./2)
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
colormap(viridis), caxis([0 1]/2), colorbar
xlim([-20 20]./2),ylim([-10 10]./2)
xlabel('{\it{l}} in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
ylabel('k in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
title(['punched: paf\_karen',sprintf('(data2, %d, %3.2f, %3.2f, Q, [%3.2f, %3.2f, %3.2f], %d)',...
                [box_size,thr_sigm,fill_flg,delta_qB,delta_sh]),newline,...
      'h = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dlq),' r.l.u.'],...
      'interpreter','latex')
axis square, grid on
set(gca,'fontsize',17,'xtick',-100:10:100,'ytick',-100:5:100,'linewidth',1,...
    'xcolor','k','ycolor','k')
% savepng([save_dirr,'plots/'],['orthoiii2_700_comparizon','_kl'],500)

%% h-l plot
n = 0; % q value along the integration axis
dlq=0.025; % half of the integration size

figure(202)
clf, dime(45,20)
subplot(1,2,1)
[dp,xp,yp,zp]=slice_cube(data1,qh,qk,ql,[hmin hmax],[-dlq dlq]+n,[lmin lmax],[0 1 0]);
imagesc(yp,xp,dp)
colormap(viridis), caxis([0 1]/2), colorbar
xlim([-20 20]/2),ylim([-10 10]/2)
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
colormap(viridis), caxis([0 1]/2), colorbar
xlim([-20 20]/2),ylim([-10 10]/2)
xlabel('{\it{l}} in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
ylabel('h in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
title(['punched: paf\_karen',sprintf('(data2, %d, %3.2f, %3.2f, Q, [%3.2f, %3.2f, %3.2f], %d)',...
                [box_size,thr_sigm,fill_flg,delta_qB,delta_sh]),newline,...
       'k = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dlq),' r.l.u.'],...
       'interpreter','latex')
axis square, grid on
set(gca,'fontsize',17,'xtick',-100:10:100,'ytick',-100:5:100,'linewidth',2,...
    'xcolor','k','ycolor','k')
% savepng([save_dirr,'plots/'],['orthoiii2_700_comparizon','_hl'],500)

%% h-k plot
n = 0; % q value along the integration axis
dlq=0.025; % half of the integration size

figure(203)
clf, dime(45,20)
subplot(1,2,1)
[dp,xp,yp,zp]=slice_cube(data1,qh,qk,ql,[hmin hmax],[kmin kmax],[-dlq dlq]+n,[0 0 1]);
imagesc(yp,xp,dp)
colormap(viridis), caxis([0 1]/2), colorbar
xlim([-5 5]), ylim([-5 5])
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
xlim([-5 5]), ylim([-5 5])
xlabel('k in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
ylabel('h in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
title(['punched: paf\_karen',sprintf('(data2, %d, %3.2f, %3.2f, Q, [%3.2f, %3.2f, %3.2f], %d)',...
                [box_size,thr_sigm,fill_flg,delta_qB,delta_sh]),newline,...
       '{\it{l}} = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dlq),' r.l.u.'],...
       'interpreter','latex')
axis square, grid on
set(gca,'fontsize',17,'xtick',-100:5:100,'ytick',-100:5:100,'linewidth',1,...
    'xcolor','k','ycolor','k')
% savepng([save_dirr,'plots/'],['orthoiii2_700_comparizon','_hk'],500)

