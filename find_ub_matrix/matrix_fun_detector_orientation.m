function f = matrix_fun_detector_orientation(x)
% function f for find_ub.m
% parameters

% x orientation matrix in real space [3 x 3]

% input - reflection list with phi, xpix, ypix, h, k, l
% note: matlab exchanges x and y on the graphics output
% keep original - so values read from matlab cursor need to be interchanged

% this means
% if the image is visualized in matlab with imagesc() function
% the 'data cursor' button used every pixel has an 'X' and 'Y'
% coordinate. The first value of the 'list_in' array is the 'Y' 
% coordinate, the second value is the 'X' coordinate. The third
% values indicates the image number (adjust the stepwidth in line 101)
% Fourth value is a dummy intensity not used followed by h, k, l, 
% where the order matters.
% The origin of the image is determined similarly as the first two 
% values of 'list_in'.
% Note: when the pixel (X,Y) is addressed directly in the array it 
% is e.g. img(Y,X).
%

load('pars_detector_orientation.mat','pars_do');
dx                 = pars_do(1);
dy                 = pars_do(2);
lambda             = pars_do(3);
ki_vec             = pars_do(4:6);
det_x              = pars_do(7:9);
det_y              = pars_do(10:12);
abc0=[pars_do(end-8:end-6);pars_do(end-5:end-3);pars_do(end-2:end)];

orgx=x(1);
orgy=x(2);
det_dist=x(3);
hor_angle=x(4);
vert_angle=x(5);

% read input list
% x-pixel y-pixel image_number h k l
filename = 'newlist.dat';
list_in=importdata(filename);

phi =list_in(:,3)*pi/180;
m=list_in(:,1);
n=list_in(:,2);

myzero=zeros(size(list_in,1),1);
myone=myzero+1;
clear('roti');
roti(1,:,:) = reshape([ cos(phi'); myzero';-sin(phi')],1,3,size(list_in,1));
roti(2,:,:) = reshape([   myzero';  myone';   myzero'],1,3,size(list_in,1));
roti(3,:,:) = reshape([ sin(phi'); myzero'; cos(phi')],1,3,size(list_in,1));

det_rot_x=[1,0,0; 0,cos(hor_angle),-sin(hor_angle); 0,sin(hor_angle),cos(hor_angle)];
det_rot_y=[cos(vert_angle),0,-sin(vert_angle); 0,1,0; sin(vert_angle),0,cos(vert_angle)];
det_rot_z=[det_x; det_y; 0,0,1];

det_matrix(:,1)=(m-orgx)*dx;
det_matrix(:,2)=(n-orgy)*dy;
det_matrix(:,3)=0;

kaa_x = det_rot_z(1,1)*det_matrix(:,1) + det_rot_z(1,2)*det_matrix(:,2) + det_rot_z(1,3)*det_matrix(:,3);
kaa_y = det_rot_z(2,1)*det_matrix(:,1) + det_rot_z(2,2)*det_matrix(:,2) + det_rot_z(2,3)*det_matrix(:,3);
kaa_z = det_rot_z(3,1)*det_matrix(:,1) + det_rot_z(3,2)*det_matrix(:,2) + det_rot_z(3,3)*det_matrix(:,3);

kaa_x_roty = det_rot_y(1,1)*kaa_x + det_rot_y(1,2)*kaa_y + det_rot_y(1,3)*kaa_z;
kaa_y_roty = det_rot_y(2,1)*kaa_x + det_rot_y(2,2)*kaa_y + det_rot_y(2,3)*kaa_z;
kaa_z_roty = det_rot_y(3,1)*kaa_x + det_rot_y(3,2)*kaa_y + det_rot_y(3,3)*kaa_z;

kaa_x_rotxy = det_rot_x(1,1)*kaa_x_roty + det_rot_x(1,2)*kaa_y_roty + det_rot_x(1,3)*kaa_z_roty;
kaa_y_rotxy = det_rot_x(2,1)*kaa_x_roty + det_rot_x(2,2)*kaa_y_roty + det_rot_x(2,3)*kaa_z_roty;
kaa_z_rotxy = det_rot_x(3,1)*kaa_x_roty + det_rot_x(3,2)*kaa_y_roty + det_rot_x(3,3)*kaa_z_roty;

kx=kaa_x_rotxy;
ky=kaa_y_rotxy;
kz=kaa_z_rotxy+det_dist;

kf_length = (sqrt(kx.^2+ky.^2+kz.^2))*lambda;
kf = permute(cat(3,(kx./kf_length)',(ky./kf_length)',(kz./kf_length)'),[1 3 2]);
ki=repmat(ki_vec,[1 1 size(list_in,1),1]);
dki=kf-ki;

hkli=permute(list_in(:,4:6),[1 2 3]);

abc = mmat(roti, repmat(abc0,[1,1,size(list_in,1)]),[1 2]);

f = sqrt(sum(sum((squeeze(abs(mmat(dki,abc,[1 2])))'-abs(hkli)).^2)))
end
