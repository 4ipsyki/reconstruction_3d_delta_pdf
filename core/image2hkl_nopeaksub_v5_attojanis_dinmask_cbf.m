function [Int_total,z_total] = image2hkl_nopeaksub_v5_attojanis_dinmask_cbf(...
                                        framenr,dirr,file,...
                                        rot_angle,...
                                        x_size,y_size,rot_ax,ub,dk,P,...
                                        mask,hmax,hmin,kmax,kmin,...
                                        lmax,lmin,n_steps,mybkg,myabs)
%                                          lmax,lmin,n_steps)
hdelta=(hmax-hmin)/(n_steps-1);
kdelta=(kmax-kmin)/(n_steps-1);
ldelta=(lmax-lmin)/(n_steps-1);

dinamic_mask1=zeros(x_size,y_size);
% dinamic_mask2=zeros(x_size,y_size);
rr1=500; % radius of the masking circle
% rr2=5; % radius of the masking circle
% rr2=1; % radius of the masking circle; If rr2=1, than only the hot pixel is masked
                                     % If rr2=2, than 3x3 square is formed
                                     % If rr2=3, than 5x5 square is formed
                                     % If rr2=4, than two 5x7 rectangles at 90 deg. are formed, resembling a circle
                                     % ...
[xx,yy]=meshgrid(1:y_size,1:x_size);


Int_total = zeros(int16(n_steps),int16(n_steps),int16(n_steps));
z_total   = zeros(int16(n_steps),int16(n_steps),int16(n_steps));

% disabling the TIFF warning
warning('off','MATLAB:imagesci:rtifc:missingPhotometricTag')

framenr=framenr(:)';
%Loop through all images
for idx=1:size(framenr,2)
    disp(framenr(idx))

%     image = double(get_pe_new5(dirr, file, framenr(idx)));
    image = permute(double(read_cbf([dirr, file, sprintf('%05d.cbf',framenr(idx))]).data),[2 1]);

    
    %     % Background on the border of each image depends both on the absorption
    %     % of the sample (which depends on the shape mostly), and the compton
    %     % background which we assume to be constant
    
    hotpix=find(image.*(~mask)>1e5);
    disp([num2str(length(hotpix)),' hot pixel(s)'])
    for hh=(hotpix(:))'
%         ix=(-2000:2000)+yy(hh); p1=find(ix>0,1,'first'); p2=find(ix<y_size+1,1,'last');
%         iy=(-10:10)+xx(hh); p3=find(iy>0,1,'first'); p4=find(iy<x_size+1,1,'last');
%         dinamic_mask1(ix(p1:p2),iy(p3:p4))=1; % masking stripes from overexposed pixels for 1 frame only
        
        circ=(xx-xx(hh)).^2+(yy-yy(hh)).^2;
        dinamic_mask1(circ<rr1^2)=1; % masking the overexposed pixels + nearby region for 1 image
%         circ=(xx-xx(hh)).^2+(yy-yy(hh)).^2;
%         dinamic_mask2(circ<rr2^2)=2000; % masking the overexposed pixels for N images
    end
%     mask_tot=mask+dinamic_mask1+dinamic_mask2; % updating the mask with the dinamic one
    try
        mask_tot=mask+dinamic_mask1; % updating the mask with the dinamic one
    catch
        disp(size(mask_tot))
        disp(size(mask))
        disp(size(dinamic_mask1))
    end
    
    mask_tot=boolean(mask_tot);
    
    image = double(image*myabs(idx)-mybkg);
    % scaling the image by the absorption factor
    % (normilized to the first image) and subtracting
    % the simulated compton scattering
    
    % phi rotation
%     phi=-(filenumber-1)*oszillation_range*pi/180+start_angle*pi/180;
    phi=-rot_angle(idx)*pi/180;
    
    
    a1=cos(phi)+rot_ax(1)^2*(1-cos(phi));
    a2=rot_ax(1)*rot_ax(2)*(1-cos(phi))-rot_ax(3)*sin(phi);
    a3=rot_ax(1)*rot_ax(3)*(1-cos(phi))+rot_ax(2)*sin(phi);
    b1=rot_ax(1)*rot_ax(2)*(1-cos(phi))+rot_ax(3)*sin(phi);
    b2=cos(phi)+rot_ax(2)^2*(1-cos(phi));
    b3=rot_ax(2)*rot_ax(3)*(1-cos(phi))-rot_ax(1)*sin(phi);
    c1=rot_ax(3)*rot_ax(1)*(1-cos(phi))-rot_ax(2)*sin(phi);
    c2=rot_ax(3)*rot_ax(2)*(1-cos(phi))+rot_ax(1)*sin(phi);
    c3=cos(phi)+rot_ax(3)^2*(1-cos(phi));
    
    rot_phi=[a1,a2,a3;b1,b2,b3;c1,c2,c3];
    
    ub_rot_phi=rot_phi*ub;
    
    % Zeile mal Spalte
    h=dk(:,:,1).*ub_rot_phi(1,1)+dk(:,:,2).*ub_rot_phi(2,1)+dk(:,:,3).*ub_rot_phi(3,1);
    k=dk(:,:,1).*ub_rot_phi(1,2)+dk(:,:,2).*ub_rot_phi(2,2)+dk(:,:,3).*ub_rot_phi(3,2);
    l=dk(:,:,1).*ub_rot_phi(1,3)+dk(:,:,2).*ub_rot_phi(2,3)+dk(:,:,3).*ub_rot_phi(3,3);
    
    pix_h = round(h/hdelta -(hmin-hdelta)/hdelta); %.*A.*B.*C;
    pix_k = round(k/kdelta -(kmin-kdelta)/kdelta); %.*A.*B.*C;
    pix_l = round(l/ldelta -(lmin-ldelta)/ldelta); %.*A.*B.*C;
    
    % apply LP-corection here
    Intens = image ./P;
    %
    final = [pix_h(1:x_size*y_size)', pix_k(1:x_size*y_size)',...
        pix_l(1:x_size*y_size)', Intens(1:x_size*y_size)',...
        mask_tot(1:x_size*y_size)'];
    
    for nn=1:size(final,1)
        if final(nn,1)<n_steps && final(nn,1)>0       && ...
                final(nn,2)>0       && final(nn,2)<n_steps && ...
                final(nn,3)<n_steps && final(nn,3)>0       && ...
                final(nn,5)==0
            Int_total(final(nn,1),final(nn,2),final(nn,3)) = ...
                Int_total(final(nn,1),final(nn,2),final(nn,3))+final(nn,4);
            z_total(final(nn,1),final(nn,2),final(nn,3))   = ...
                z_total(final(nn,1),final(nn,2),final(nn,3))+1;
        end
    end
    dinamic_mask1=dinamic_mask1*0; % reinitializing to zero
%     idxm=~~dinamic_mask2; % having 1 if masked and 0 if not. Maybe no need for this.
%     dinamic_mask2(idxm)=dinamic_mask2(idxm)-1; % decreasing the mask value until it reaches 0 and thus it won't be masked anymore
end

% restoring the TIFF warning
warning('on','MATLAB:imagesci:rtifc:missingPhotometricTag')

end
