base_dirr        = '/asap3/petra3/gpfs/p07/2019/data/11006746/';
name_ip          = 'ip_700_allpos_02';
pars_dirr        = [base_dirr,'processed/inputpar/'];

addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/core/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/additional/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/find_ub_matrix/']);
load([pars_dirr,name_ip,'.mat']);

% saving directory
if ~exist(save_dirr, 'dir')
       mkdir(save_dirr)
end

clear mask; load([save_dirr,'masks_700.mat']);
clear abs_corr; load([save_dirr,'abs_corr_700.mat']);

% the following loop is made to reconstruct and save (in a temporary folder) partial
% files corresponding to the six different lateral positions of the detector. each of
% these six (macro)positions is composed of four (micro)positions to cover the gaps of
% the detector. Use merge.m to merge the reconstructions (after inspecting them separately).
% modify this part accordingly to your data-acquisition scheme.
% For a one combined (several runs) reconstruction comment the first for loop, and
% uncomment the relevant commented parts (and commenting the respective lines just below)

for macropos=1:6
    Int_total = zeros(n_steps);
    z_total = zeros(n_steps);
    % for nr=1:size(dirr,2)
    for minipos=1:4
        nr=4*(macropos-1)+minipos;
        % myab=ones(1,size(frame_nr{nr},2)); % use this if no absorption correction is needed
        myab=1./abs_corr{nr}; % absorption and/or monitor corrections. Add here if aplicable.
        % reconstruction
        [dk,P]=calc_dk_P(orgx(nr),orgy(nr),det_dist(nr),ki_vec,rot_ax,...
                         pol_plane_normal,pol_deg,x_size,y_size,...
                         hor_angle(nr),det_x,det_y,vert_angle(nr),dx,dy,lambda);

        % [Int_total,z_total]=image2hkl_v6(...
        [Int_partial,z_partial]=image2hkl_v6(...
                                frame_nr{nr},...
                                dirr{nr},file{nr},rot_angle{nr},...
                                x_size,y_size,rot_ax,ub,dk,P,boolean(mask{nr}),...
                                hmax,hmin,kmax,kmin,lmax,lmin,n_steps,0,myab,...
                                data_type,din_mask_pars);
        Int_total = Int_total + Int_partial;
        z_total = z_total + z_partial;
    end
    % if the reconstruction is not partial, the normilized intensities can be directly saved
    % Int_total_t = Int_total./z_total;
    % Int_total_t(isnan(Int_total_t)) = 0;
    % saving
    file_name=[strrep(file{nr}(1:end-1),'_pos4',''),'_dinmask_part',num2str(macropos),'.mat'];
    % save([save_dirr,file_name],'Int_total_t',...
    save([save_dirr,'rec_700mm_temp/',file_name],'Int_total','z_total',...
                                                   'hmin','hmax',...
                                                   'kmin','kmax',...
                                                   'lmin','lmax',...
                                                   'ub','-v7.3');
end

close all
clear all
exit()
