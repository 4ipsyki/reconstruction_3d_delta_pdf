base_dirr        = '/asap3/petra3/gpfs/p07/2019/data/11006746/';
name_ip          = 'ip_700_allpos_02';
pars_dirr        = [base_dirr,'processed/inputpar/'];

addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/core/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/additional/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/find_ub_matrix/']);
load([pars_dirr,name_ip,'.mat']);

loadnm=[
'orthoiii2_00007_dinmask_part6.mat';
'orthoiii2_00006_dinmask_part5.mat';
'orthoiii2_00005_dinmask_part4.mat';
'orthoiii2_00003_dinmask_part3.mat';
'orthoiii2_00002_dinmask_part2.mat';
'orthoiii2_00001_dinmask_part1.mat'];

I_part=zeros(n_steps);
z_part=zeros(n_steps);

for ii=1:size(loadnm,1)
    save_name = [save_dirr,'rec_700mm_temp_2/',loadnm(ii,:)];
    load(save_name,'Int_total','z_total','hmin','hmax','kmin','kmax','lmin','lmax');

    I_part=I_part+Int_total;
    z_part=z_part+z_total;
end
Int_total_t=I_part./z_part;
Int_total_t(isnan(Int_total_t))=0;
save_name = [save_dirr,'orthoiii2_700_merged_v2.mat'];
% rescaling axes since from plotting the peaks are more aligned with the
% grid. The sample is very twinned in h-k planes
hmin=hmin/1.007;hmax=hmax/1.007;
kmin=kmin*1.009;kmax=kmax*1.009;
ub_scaled=[ub(:,1)*1.007, ub(:,2)/1.009, ub(:,3)];
save(save_name,'Int_total_t','hmin','hmax','kmin','kmax','lmin','lmax',...
    'ub','ub_scaled','-v7.3');

