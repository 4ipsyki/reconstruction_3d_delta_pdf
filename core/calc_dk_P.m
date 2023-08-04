function [dk,P] = calc_dk_P_v1(orgx,orgy,det_dist,ki_vec,rot_ax,...
                            pol_plane_normal,pol_deg,x_size,y_size,...
                            hor_angle,det_x,det_y,vert_angle,dx,dy,lambda)
% general calcs
m=repmat((1:x_size)',1,y_size);
n=repmat( 1:y_size  ,x_size,1);

det_rot_x=[1,0,0; 0,cos(hor_angle),-sin(hor_angle); 0,sin(hor_angle),cos(hor_angle)];     % rotation around horizontal axis
det_rot_y=[cos(vert_angle),0,-sin(vert_angle); 0,1,0; sin(vert_angle),0,cos(vert_angle)]; % rotation around vertical axis
det_rot_z=[det_x; det_y; 0,0,1];

det_matrix(:,:,1)=(m-orgx)*dx;
det_matrix(:,:,2)=(n-orgy)*dy;
det_matrix(:,:,3)=zeros(x_size,y_size);

kaa_x = det_rot_z(1,1)*det_matrix(:,:,1) + det_rot_z(1,2)*det_matrix(:,:,2) + det_rot_z(1,3)*det_matrix(:,:,3);
kaa_y = det_rot_z(2,1)*det_matrix(:,:,1) + det_rot_z(2,2)*det_matrix(:,:,2) + det_rot_z(2,3)*det_matrix(:,:,3);
kaa_z = det_rot_z(3,1)*det_matrix(:,:,1) + det_rot_z(3,2)*det_matrix(:,:,2) + det_rot_z(3,3)*det_matrix(:,:,3);

kaa_x_roty = det_rot_y(1,1)*kaa_x + det_rot_y(1,2)*kaa_y + det_rot_y(1,3)*kaa_z;
kaa_y_roty = det_rot_y(2,1)*kaa_x + det_rot_y(2,2)*kaa_y + det_rot_y(2,3)*kaa_z;
kaa_z_roty = det_rot_y(3,1)*kaa_x + det_rot_y(3,2)*kaa_y + det_rot_y(3,3)*kaa_z;

kaa_x_rotxy = det_rot_x(1,1)*kaa_x_roty + det_rot_x(1,2)*kaa_y_roty + det_rot_x(1,3)*kaa_z_roty;
kaa_y_rotxy = det_rot_x(2,1)*kaa_x_roty + det_rot_x(2,2)*kaa_y_roty + det_rot_x(2,3)*kaa_z_roty;
kaa_z_rotxy = det_rot_x(3,1)*kaa_x_roty + det_rot_x(3,2)*kaa_y_roty + det_rot_x(3,3)*kaa_z_roty;

kx=kaa_x_rotxy;
ky=kaa_y_rotxy;
kz=kaa_z_rotxy+det_dist;

kf=cat(3,kx,ky,kz).*repmat(1./(sqrt(kx.^2+ky.^2+kz.^2))*1/lambda,[1 1 3]);
ki=cat(3,repmat(ki_vec(1),x_size,y_size),repmat(ki_vec(2),x_size,y_size),repmat(ki_vec(3),x_size,y_size));
dk=kf-ki;

% LP-correction
% Kabsch J. Appl. Cryst. (1988) 21, 916-924:
% Kabsch Acta Crystallogr D Biol Crystallogr. 2014, 70 2204–2216: 
% L=|u dot (S x S0)| / (|S| dot |S0|)
% L correction is not performed here since it is included when the summed intensity in each voxel is normilized to the number of entries: Int_total./z_total

% P=(1-2p)[1-(n dot S / |S|)^2] + p{1+[S dot S0 / (|S||S0|)]^2}
% S0: incident beam
% S: diffracted beam
% u: unit vector along the direction of the rotation axis
% p: degree of polarization, 0.5=unpolarized
% n: polarization plane normal
norm_ki=sqrt(ki(:,:,1).^2+ki(:,:,2).^2+ki(:,:,3).^2);
norm_kf=sqrt(kf(:,:,1).^2+kf(:,:,2).^2+kf(:,:,3).^2);

cross_kf_ki(:,:,1) = kf(:,:,2).*ki(:,:,3) - kf(:,:,3).*ki(:,:,2);
cross_kf_ki(:,:,2) = kf(:,:,3).*ki(:,:,1) - kf(:,:,1).*ki(:,:,3);
cross_kf_ki(:,:,3) = kf(:,:,1).*ki(:,:,2) - kf(:,:,2).*ki(:,:,1);

dot_rot_ax_cross_kf_ki=rot_ax(1).*cross_kf_ki(:,:,1) + rot_ax(2).*cross_kf_ki(:,:,2) + rot_ax(3).*cross_kf_ki(:,:,3);

dot_pol_plane_normal_kf = pol_plane_normal(1)*kf(:,:,1) + pol_plane_normal(2)*kf(:,:,2) + pol_plane_normal(3)*kf(:,:,3);
dot_kf_ki = kf(:,:,1).*ki(:,:,1) + kf(:,:,2).*ki(:,:,2) + kf(:,:,3).*ki(:,:,3); 

P=(1-2*pol_deg).*(1-(dot_pol_plane_normal_kf ./ norm_kf).^2) + pol_deg*(1+(dot_kf_ki ./ (norm_kf .* norm_ki)).^2);

end

