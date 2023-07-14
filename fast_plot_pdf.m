base_dirr        = '/asap3/petra3/gpfs/p07/2019/data/11006746/';
name_ip          = 'ip_700_allpos_02';
pars_dirr        = [base_dirr,'processed/inputpar/'];

addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/core/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/additional/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/find_ub_matrix/']);
load([pars_dirr,name_ip,'.mat'],'save_dirr');

loadnm=['orthoiii2_700_merged_v2_punched_box92_sg3_bragg0p27_0p6_0p6_filled_sg0p3_kernel5_sym_ifft_elipse.mat'];
load([save_dirr,loadnm],'Int_total_t','dx','dy','dz');

Int_total_t(isnan(Int_total_t))=0;
data1=Int_total_t;

%% k-l plot
n = 0; % q value along the integration axis
dld=0.05; % half of the integration size

figure(301)
clf, dime(20,20)
[dp,xp,yp,zp]=slice_cube(data1,dx,dy,dz,[-dld dld]+n,[min(dy) max(dy)],[min(dz) max(dz)],[1 0 0]);
imagesc(yp,xp,dp), view([0 -90])
colormap(redblue), caxis([-1 1]*1e-5), colorbar
xlim([-0.7 5.7]), ylim([-0.7 5.7])
xlabel('{\it{w}} in ({\it{u}},{\it{v}},{\it{w}}) [l.u.]','interpreter','latex')
ylabel('{\it{v}} in ({\it{u}},{\it{v}},{\it{w}}) [l.u.]','interpreter','latex')
title([strrep(loadnm(1:60),'_','\_'),newline,strrep(loadnm(61:end),'_','\_'),newline,...
      '{\it{u}} = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dld),' l.u.'],...
      'interpreter','latex')
axis square, grid on
set(gca,'fontsize',15,'xtick',-20:1:20,'ytick',-20:1:20,'linewidth',1,...
    'xcolor','k','ycolor','k')
% savepng([save_dirr,'plots/'],[loadnm,'_vw'],500)

%% h-l plot
n = 0; % q value along the integration axis
dld=0.05; % half of the integration size

figure(302)
clf, dime(20,20)
[dp,xp,yp,zp]=slice_cube(data1,dx,dy,dz,[min(dx) max(dx)],[-dld dld]+n,[min(dz) max(dz)],[0 1 0]);
imagesc(yp,xp,dp), view([0 -90])
colormap(redblue), caxis([-1 1]*1e-5), colorbar
xlim([-0.7 5.7]), ylim([-0.7 5.7])
xlabel('{\it{w}} in ({\it{u}},{\it{v}},{\it{w}}) [l.u.]','interpreter','latex')
ylabel('{\it{u}} in ({\it{u}},{\it{v}},{\it{w}}) [l.u.]','interpreter','latex')
title([strrep(loadnm(1:60),'_','\_'),newline,strrep(loadnm(61:end),'_','\_'),newline,...
      '{\it{v}} = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dld),' l.u.'],...
      'interpreter','latex')
axis square, grid on
set(gca,'fontsize',15,'xtick',-20:1:20,'ytick',-20:1:20,'linewidth',1,...
    'xcolor','k','ycolor','k')
% savepng([save_dirr,'plots/'],[loadnm,'_uw'],500)

%% h-k plot
n = 0; % q value along the integration axis
dld=0.05; % half of the integration size

figure(303)
clf, dime(20,20)
[dp,xp,yp,zp]=slice_cube(data1,dx,dy,dz,[min(dx) max(dx)],[min(dy) max(dy)],[-dld dld]+n,[0 0 1]);
imagesc(yp,xp,dp), view([0 -90])
colormap(redblue), caxis([-1 1]*1e-5), colorbar
xlim([-0.7 5.7]), ylim([-0.7 5.7])
xlabel('{\it{v}} in ({\it{u}},{\it{v}},{\it{w}}) [l.u.]','interpreter','latex')
ylabel('{\it{u}} in ({\it{u}},{\it{v}},{\it{w}}) [l.u.]','interpreter','latex')
title([strrep(loadnm(1:60),'_','\_'),newline,strrep(loadnm(61:end),'_','\_'),newline,...
      '{\it{w}} = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dld),' l.u.'],...
      'interpreter','latex')
axis square, grid on
set(gca,'fontsize',15,'xtick',-20:1:20,'ytick',-20:1:20,'linewidth',1,...
    'xcolor','k','ycolor','k')
% savepng([save_dirr,'plots/'],[loadnm(1:end-4),'_uv'],500)

