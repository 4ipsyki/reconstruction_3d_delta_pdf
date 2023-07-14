base_dirr        = '/asap3/petra3/gpfs/p07/2019/data/11006746/';
name_ip          = 'ip_700_allpos_02';
pars_dirr        = [base_dirr,'processed/inputpar/'];

addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/core/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/additional/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/find_ub_matrix/']);
load([pars_dirr,name_ip,'.mat'],'save_dirr');

myfiles={'orthoiii2_700_merged_v2_punched_box92_sg3_bragg0p27_0p6_0p6_filled_sg0p3_kernel5_sym.mat';};

idxfile=1;
loadnm=myfiles{idxfile};
load([save_dirr,loadnm],'Int_total_t','hmin','hmax','kmin','kmax','lmin','lmax','ub');

data1=Int_total_t; data1(isnan(data1))=0; clear Int_total_t
n_steps=size(data1);
qh=linspace(hmin,hmax,n_steps(1));
qk=linspace(kmin,kmax,n_steps(2));
ql=linspace(lmin,lmax,n_steps(3));

%% elipsoid/circle high-Q cutoff
[q2,q1,q3]=meshgrid(qk,qh,ql);

alp=norm(ub(:,1)); blp=norm(ub(:,2)); clp=norm(ub(:,3));

% maximum Q values (in r.l.u.) for each direction
qhmax=15.5;
qkmax=13.5;
qlmax=40.1;

% % making it circular in the hkl space
% Qmax = [    qhmax       qhmax/hmax*kmax   qhmax/hmax*lmax;
%         qkmax/kmax*hmax       qkmax       qkmax/kmax*lmax;
%         qlmax/lmax*hmax   qlmax/lmax*kmax       qlmax];
% [~,idx_min]=min(Qmax(:,1));
% qhmax=Qmax(idx_min,1);
% qkmax=Qmax(idx_min,2);
% qlmax=Qmax(idx_min,3);

idx_highq=(q1.^2/qhmax^2 + q2.^2/qkmax^2 + q3.^2/qlmax^2) > 1;

% making last high-Q data going to zero intensity by substructing a median
% intensity determined from the border of the data
idx_1=(q1.^2/(qhmax-0.1)^2 + q2.^2/(qkmax-0.1)^2 + q3.^2/(qlmax-0.1)^2) > 1;
idx_b=boolean(~idx_highq.*idx_1);
int_b=data1(idx_b); int_b(int_b==0)=[]; int_b(isnan(int_b))=[];
bck=median(int_b);
data1=data1-bck;
% resetting no-data region to zero
data1(idx_highq)=0;

%% iFFT
data1(isnan(data1))=0;
data2=ifftshift(real(ifftn(fftshift(data1))));
dx = linspace(-n_steps(1)/2,n_steps(1)/2,n_steps(1))*(1/(2*hmax));
dy = linspace(-n_steps(2)/2,n_steps(2)/2,n_steps(2))*(1/(2*kmax));
dz = linspace(-n_steps(3)/2,n_steps(3)/2,n_steps(3))*(1/(2*lmax));

%% saving
Int_total_t=data2;
Int_total_t(isnan(Int_total_t))=0;
save_name=[loadnm(1:end-4),'_ifft_elipse.mat'];
% save_name=[loadnm(1:end-4),'_ifft_circle.mat'];
save([save_dirr,save_name],'Int_total_t','dx','dy','dz','ub','-v7.3');

