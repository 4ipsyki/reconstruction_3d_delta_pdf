% name of the input par file to be saved
name_ip          = 'ip_700_allpos_02';

% main directory
base_dirr        = '/asap3/petra3/gpfs/p07/2019/data/11006746/';
% saving directory for reconstructions
save_dirr        = [base_dirr,'processed/reconstructions/ybco_oleh/'];
% saving directory for input-prameter file (containing all these pars)
pars_dirr        = [base_dirr,'processed/inputpar/'];

% sub directory for data
main_dirr        = 'raw/pilatus/';
% generic file name for frames
main_name        = 'orthoiii2_pos';

% detector type and orientation:
% 1 - Pilatus2M;                2 - PerkinElmer @P21.1;
% 3 - PerkinElmer @P07/EH2      4 - Custom --> fill parameters section
detector         = 1;
% extension of the frames
data_type        = 1; % 1 - TIF; 2 - CBF;

% scan number(s) the different runs to be reconstrcuted
scan_num         = [1:3,5:7];

% scaling factor for angle reading (if necessary)
angle_conv       = 1;
% logging of angle for each frame: 0 - NO; 1 - YES;
% if YES, add/modify log-reading script (se below)
angle_logs       = 0;
% angle step per frame in degrees (not used if angle_logs = 1)
angle_step       = 0.02;
% total number of frames per scan (not used if angle_logs = 1)
tot_frames       = 10002;

% data path / filename / angle-per-frame / frame-number custom-code generator
% modify this part according to your data
clear dirr file rot_angle framenr
ll = length(scan_num);
dirr = cell(ll,1); file = cell(ll,1);
rot_angle = cell(ll,1); frame_nr = cell(ll,1);
idx = 1;
for nr = scan_num
    for mypos = 1:4
        name = [main_name,sprintf('%d_%05d',[mypos,nr])];
        dirr{idx} = [base_dirr,main_dirr,name,'/'];
        file{idx} = [name,'_'];
        if angle_logs==1
            nm_attojanis = [base_dirr,main_dirr,'log/',name,'.dat'];
            [rot_angle{idx},frame_nr{idx}] = read_log_attojanis_v2(nm_attojanis);
            rot_angle{idx} = rot_angle{idx}*angle_conv;
        else
            rot_angle{idx} = (0:1:(tot_frames-1))*angle_step*angle_conv;
            frame_nr{idx}  = 1:1:tot_frames;
        end
        idx = idx+1;
    end
end
clear scan_num mypos ll idx name nr nm_attojanis main_name main_dirr tot_frames angle_conv angle_step angle_logs

% for generation of ub matrix
% number of frames to be considered for ub matrix
n_pics           = 10000;
if n_pics>length(frame_nr{1}); n_pics = length(frame_nr{1}); end
% nominal lattice parameters
lattice          = [3.9 3.9 11.94 90 90 90];

% 10000 frames and 700mm data
% f =    3.6888    3.8468    3.8468   11.7071   90.0000   90.0000   90.0000
ub                 = [  0.198664081781044,  0.277196592085856, 11.661013808659812;
                        3.840046348727848,  0.098616305270907, -0.627623679330360;
                       -0.113088647012164,  3.835576787592934, -0.826603020939706];
% reconstruction properties
% number of voxels for all three directions
n_steps          = [920 920 920]; if length(n_steps)~=3; n_steps=ones(1,3)*n_steps(1); end
% limits in each direction
hmin             = -18.38; % res. of 0.04 rlu
hmax             =  18.38;
kmin             = -18.38;
kmax             =  18.38;
lmin             = -45.95; % res. of 0.1 rlu
lmax             =  45.95;
% --> resolution = (q_max-q_min)/(n_steps-1);

% dinamic masking parameters [threshold
%                             radius_of_blooming_masking_circle_px
%                             radius_of_hot_pixel_masking_circle_px
%                             masking_frames_number_for_hot_pixel]
% if threshold==0  no dinamic masking is executed
% masking_frames_number_for_hot_pixel is the number of successive frames on which
% radius_of_hot_pixel_masking_circle_px will be applied. This is necessary in case
% of detectors with after-glow problem
din_mask_pars    = [1e4 300 2 1];

% det_dist, orgx, orgy, hor_angle and vert_angle are obtained using:
% poni2org([[base_dirr,'processed/pilatus/ce02_8/','ceo2_8_pos1_00001.poni'];...
%           [base_dirr,'processed/pilatus/ce02_8/','ceo2_8_pos2_00001.poni'];...
%           [base_dirr,'processed/pilatus/ce02_8/','ceo2_8_pos3_00001.poni'];...
%           [base_dirr,'processed/pilatus/ce02_8/','ceo2_8_pos4_00001.poni'];...
%           [base_dirr,'processed/pilatus/ce02_8/','ceo2_8_pos1_00002.poni'];...
%           [base_dirr,'processed/pilatus/ce02_8/','ceo2_8_pos2_00002.poni'];...
%           [base_dirr,'processed/pilatus/ce02_8/','ceo2_8_pos3_00002.poni'];...
%           [base_dirr,'processed/pilatus/ce02_8/','ceo2_8_pos4_00002.poni']],dx);
% missing calibrations are filled by extrapolation from the existing ones
det_dist   = [ 647.7194  647.7695  647.6945  647.6373...
               644.7823  644.8045  644.7327  644.7120...
               641.8452  641.8395  641.7709  641.7867...
               647.7194  647.7695  647.6945  647.6373...
               644.7823  644.8045  644.7327  644.7120...
               641.8452  641.8395  641.7709  641.7867];   % mm
orgx       = [ 784.9402  814.0199  813.9826  784.9098...
               783.2408  812.3511  812.3214  783.2224...
               781.5414  810.6823  810.6602  781.5350...
              2180.2890 2209.3687 2209.3314 2180.2586...
              2178.5896 2207.6999 2207.6702 2178.5712...
              2176.8902 2206.0311 2206.0090 2176.8838];   % px
orgy       = [1585.1836 1585.2154 1556.1269 1556.1755...
               420.6978  420.5714  391.4581  391.5398...
              -743.7880 -744.0726 -773.2107 -773.0959...
              1585.1836 1585.2154 1556.1269 1556.1755...
               420.6978  420.5714  391.4581  391.5398...
              -743.7880 -744.0726 -773.2107 -773.0959];   % px
hor_angle  = [0.00305359 0.00313093 0.00305506 0.00299179...
              0.00249166 0.00237141 0.00225678 0.00234223...
              0.00192973 0.00161189 0.00145850 0.00169267...
              0.00305359 0.00313093 0.00305506 0.00299179...
              0.00249166 0.00237141 0.00225678 0.00234223...
              0.00192973 0.00161189 0.00145850 0.00169267];   % rad
vert_angle = [0.01396220 0.01416580 0.01404766 0.01395166...
              0.01372559 0.01377472 0.01377769 0.01386753...
              0.01348898 0.01338364 0.01350772 0.01378340...
              0.01396220 0.01416580 0.01404766 0.01395166...
              0.01372559 0.01377472 0.01377769 0.01386753...
              0.01348898 0.01338364 0.01350772 0.01378340];   % rad

% x-ray wavelength in Angstroem
lambda           = 12.3984/98.19;
% incident beam direction
ki_vec           = [0.0  0.0  1/lambda];

% rotation axis
rot_ax           = [0.0 1.0 0.0];
% beam polarization
pol_deg          = 0.920;
% polarization plane normal
pol_plane_normal = [0.000 1.000  0.000];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch detector
    case 1 % Pilatus2M
        det_x    = [ 0.00000   1.00000   0.00000];
        det_y    = [ 1.00000   0.00000   0.00000];
        x_size   = 1679; % number of pixels
        y_size   = 1475;
        dx       = 0.172; % pixel size
        dy       = 0.172;
    case 2 % PerkinElmer @P21.1
        det_x    = [ 0.00000   1.00000   0.00000]; 
        det_y    = [ 1.00000   0.00000   0.00000];
        x_size   = 2048; % number of pixels
        y_size   = 2048; 
        dx       = 0.2; % pixel size
        dy       = 0.2;
    case 3 % Perkin Elmer @P07/EH2
        det_x    = [ 1.00000   0.00000   0.00000];
        det_y    = [ 0.00000  -1.00000   0.00000];
        x_size   = 2048; % number of pixels
        y_size   = 2048; 
        dx       = 0.2; % pixel size
        dy       = 0.2;
    case 4 % Custom
        det_x    = [ 1.00000   0.00000   0.00000];
        det_y    = [ 0.00000   1.00000   0.00000];
        x_size   = 1234; % number of pixels
        y_size   = 1234; 
        dx       = 0.1; % pixel size
        dy       = 0.1;
end
clear detector

myfolder = save_dirr;
if ~exist(myfolder, 'dir')
      mkdir(myfolder)
end
clear myfolder

myfolder = pars_dirr;
if ~exist(myfolder, 'dir')
      mkdir(myfolder)
end
clear myfolder

save([pars_dirr,name_ip,'.mat'],...
                        'base_dirr','save_dirr','pars_dirr',...
                        'dirr','file','data_type',...
                        'det_dist','orgx','orgy','hor_angle','vert_angle',...                        
                        'x_size','y_size','dx','dy','lambda','det_x','det_y',...
                        'ki_vec','frame_nr','rot_angle',...
                        'rot_ax','pol_deg','pol_plane_normal',...
                        'n_steps','hmin','hmax','kmin','kmax','lmin','lmax',...
                        'din_mask_pars','n_pics','lattice','ub');
                    
clear                   dirr file data_type...
                        det_dist orgx orgy hor_angle vert_angle ...
                        x_size y_size dx dy lambda det_x det_y ...
                        ki_vec frame_nr rot_angle ...
                        rot_ax pol_deg pol_plane_normal ...
                        n_steps hmin hmax kmin kmax lmin lmax ...
                        din_mask_pars n_pics lattice ub ...

