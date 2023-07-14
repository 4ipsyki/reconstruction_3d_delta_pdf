function [ h,k,l ] = hkl_from_abc( m,n,rot_angle, abc,pars )

orgx = pars(1);
orgy = pars(2);
dx = pars(3);
dy = pars(4);
det_dist = pars(5);
lambda = pars(6);
ki_vec = [pars(7),pars(8),pars(9)];
det_x = [pars(10),pars(11),pars(12)];
det_y = [pars(13),pars(14),pars(15)];
hor_angle = pars(16);
vert_angle = pars(17);

phi=rot_angle(:)*pi/180;

myzero=0;
myone=1;
clear('roti');
roti_1 = [ cos(phi'); myzero';-sin(phi')];
roti_2 = [   myzero';  myone';   myzero'];
roti_3 = [ sin(phi'); myzero'; cos(phi')];

det_rot_x=[1,0,0; 0,cos(hor_angle),-sin(hor_angle); 0,sin(hor_angle),cos(hor_angle)];     % rotation around horizontal axis
det_rot_y=[cos(vert_angle),0,-sin(vert_angle); 0,1,0; sin(vert_angle),0,cos(vert_angle)]; % rotation around vertical axis
det_rot_z=[det_x; det_y; 0,0,1]; % 

det_matrix(1)=(m-orgx)*dx;
det_matrix(2)=(n-orgy)*dy;
det_matrix(3)=0;            % this is for rotation around beam position on detector (see below)

kaa_x = det_rot_z(1,1)*det_matrix(1) + det_rot_z(1,2)*det_matrix(2) + det_rot_z(1,3)*det_matrix(3);
kaa_y = det_rot_z(2,1)*det_matrix(1) + det_rot_z(2,2)*det_matrix(2) + det_rot_z(2,3)*det_matrix(3);
kaa_z = det_rot_z(3,1)*det_matrix(1) + det_rot_z(3,2)*det_matrix(2) + det_rot_z(3,3)*det_matrix(3);

kaa_x_roty = det_rot_y(1,1)*kaa_x + det_rot_y(1,2)*kaa_y + det_rot_y(1,3)*kaa_z;
kaa_y_roty = det_rot_y(2,1)*kaa_x + det_rot_y(2,2)*kaa_y + det_rot_y(2,3)*kaa_z;
kaa_z_roty = det_rot_y(3,1)*kaa_x + det_rot_y(3,2)*kaa_y + det_rot_y(3,3)*kaa_z;

kaa_x_rotxy = det_rot_x(1,1)*kaa_x_roty + det_rot_x(1,2)*kaa_y_roty + det_rot_x(1,3)*kaa_z_roty;
kaa_y_rotxy = det_rot_x(2,1)*kaa_x_roty + det_rot_x(2,2)*kaa_y_roty + det_rot_x(2,3)*kaa_z_roty;
kaa_z_rotxy = det_rot_x(3,1)*kaa_x_roty + det_rot_x(3,2)*kaa_y_roty + det_rot_x(3,3)*kaa_z_roty;

kx=kaa_x_rotxy;
ky=kaa_y_rotxy;
kz=kaa_z_rotxy+det_dist;   % this is for rotation around beam position on detector (see above)

kf_length = (sqrt(kx.^2+ky.^2+kz.^2))*lambda;
kf = [(kx./kf_length)',(ky./kf_length)',(kz./kf_length)'];
ki=ki_vec;
dk=kf-ki;

rot_phi=[roti_1';roti_2';roti_3'];

ub=abc;
ub_rot_phi=rot_phi*ub;
    
h=dk(1).*ub_rot_phi(1,1)+dk(2).*ub_rot_phi(2,1)+dk(3).*ub_rot_phi(3,1);
k=dk(1).*ub_rot_phi(1,2)+dk(2).*ub_rot_phi(2,2)+dk(3).*ub_rot_phi(3,2);
l=dk(1).*ub_rot_phi(1,3)+dk(2).*ub_rot_phi(2,3)+dk(3).*ub_rot_phi(3,3);

end
