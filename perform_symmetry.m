base_dirr        = '/asap3/petra3/gpfs/p07/2019/data/11006746/';
name_ip          = 'ip_700_allpos_02';
pars_dirr        = [base_dirr,'processed/inputpar/'];

addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/core/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/additional/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/find_ub_matrix/']);
load([pars_dirr,name_ip,'.mat'],'save_dirr');

myfiles={'orthoiii2_700_merged_v2.mat';};

% Laue group symmetry operations
idx_sym = [ 1  2  3; -1 -2  3; -1  2 -3; 1 -2 -3;...
           -1 -2 -3; -1  2  3;  1 -2  3; 1  2 -3];

for idxfile=1:size(myfiles,1)
    loadnm=myfiles{idxfile};
    % load([save_dirr,loadnm],'Int_total','z_total','idx_paf','hmin','hmax','kmin','kmax','lmin','lmax','ub'); 
% Int_total_t=Int_total./z_total; Int_total_t(isnan(Int_total_t))=0; Int_total_t(idx_paf)=0;
    load([save_dirr,loadnm],'Int_total','z_total','hmin','hmax','kmin','kmax','lmin','lmax','ub');
Int_total_t=Int_total./z_total; Int_total_t(isnan(Int_total_t))=0;

    data1=Int_total_t; data1(isnan(data1))=0; clear Int_total_t
    n_steps=size(data1);
    qh=linspace(hmin,hmax,n_steps(1));
    qk=linspace(kmin,kmax,n_steps(2));
    ql=linspace(lmin,lmax,n_steps(3));
    
    % applying symmetries
    Int_total(isnan(Int_total))=0; data21=apply_symmetry(Int_total,idx_sym);
    z_total(isnan(z_total))=0; data22=apply_symmetry(z_total,idx_sym);
    data2=data21./data22; data2(isnan(data2))=0;
    
    % saving
    Int_total=data21;
    z_total=data22;
    save_name=[loadnm(1:end-4),'_sym.mat'];
    % save([save_dirr,save_name],'Int_total','z_total','idx_paf','idx_sym',...
    save([save_dirr,save_name],'Int_total','z_total','idx_sym',...
                    'hmin','hmax','kmin','kmax','lmin','lmax','ub','-v7.3')
end
%% plotting
%% k-l plot
n = 12; % q value along the integration axis
dlq=0.025; % half of the integration size

% figure(101)
figure(201)
clf, dime(45,20)
subplot(1,2,1)
[dp,xp,yp,zp]=slice_cube(data1,qh,qk,ql,[-dlq dlq]+n,[kmin kmax],[lmin lmax],[1 0 0]);
imagesc(yp,xp,dp)
colormap(viridis), caxis([0 1]*0.7), colorbar
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
xlabel('{\it{l}} in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
ylabel('k in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
title(['symmetriezed',newline,...
      'h = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dlq),' r.l.u.'],...
      'interpreter','latex')
axis square, grid on
set(gca,'fontsize',17,'xtick',-100:10:100,'ytick',-100:5:100,'linewidth',1,...
    'xcolor','k','ycolor','k')
% savepng([save_dirr,'plots/'],[save_name(1:end-4),'_kl'],500)

%% h-l plot
n = 10; % q value along the integration axis
dlq=0.025; % half of the integration size

figure(202)
clf, dime(45,20)
subplot(1,2,1)
[dp,xp,yp,zp]=slice_cube(data1,qh,qk,ql,[hmin hmax],[-dlq dlq]+n,[lmin lmax],[0 1 0]);
imagesc(yp,xp,dp)
colormap(viridis), caxis([0 1]*0.7), colorbar
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
xlabel('{\it{l}} in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
ylabel('h in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
title(['symmetriezed',newline,...
       'k = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dlq),' r.l.u.'],...
       'interpreter','latex')
axis square, grid on
set(gca,'fontsize',17,'xtick',-100:10:100,'ytick',-100:5:100,'linewidth',2,...
    'xcolor','k','ycolor','k')
% savepng([save_dirr,'plots/'],[save_name(1:end-4),'_hl'],500)

%% h-k plot
n = 10; % q value along the integration axis
dlq=0.025; % half of the integration size

figure(203)
clf, dime(45,20)
subplot(1,2,1)
[dp,xp,yp,zp]=slice_cube(data1,qh,qk,ql,[hmin hmax],[kmin kmax],[-dlq dlq]+n,[0 0 1]);
imagesc(yp,xp,dp)
colormap(viridis), caxis([0 1]/2), colorbar
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
xlabel('k in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
ylabel('h in (h,k,{\it{l}}) [r.l.u.]','interpreter','latex')
title(['symmetriezed',newline,...
       '{\it{l}} = ',sprintf('%4.2f',zp),'$\pm$',sprintf('%4.3f',dlq),' r.l.u.'],...
       'interpreter','latex')
axis square, grid on
set(gca,'fontsize',17,'xtick',-100:5:100,'ytick',-100:5:100,'linewidth',1,...
    'xcolor','k','ycolor','k')
% savepng([save_dirr,'plots/'],[save_name(1:end-4),'_hk'],500)

