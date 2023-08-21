function [Int_total,z_total] = image2hkl(frame_nr,dirr,file,...
                                         rot_angle,...
                                         x_size,y_size,rot_ax,ub,dk,P,...
                                         mask,hmax,hmin,kmax,kmin,...
                                         lmax,lmin,n_steps,mybkg,myabs,...
                                         data_type,din_mask_pars)

hdelta=(hmax-hmin)/(n_steps(1)-1);
kdelta=(kmax-kmin)/(n_steps(2)-1);
ldelta=(lmax-lmin)/(n_steps(3)-1);

din_thr = din_mask_pars(1);
if din_thr > 0
    dinamic_mask1=zeros(x_size,y_size);
    dinamic_mask2=zeros(x_size,y_size);
    din_thr = din_mask_pars(1);
    din_r1  = din_mask_pars(2); % radius of the masking circle
    din_r2  = din_mask_pars(3); % radius of the masking circle
    din_fr  = din_mask_pars(4);
end

[xx,yy]=meshgrid(1:y_size,1:x_size);


Int_total = zeros(n_steps);
z_total   = zeros(n_steps);

% disabling the TIFF warning
warning('off','MATLAB:imagesci:rtifc:missingPhotometricTag')

frame_nr=frame_nr(:)';
%Loop through all images
for idx=1:size(frame_nr,1)
    disp(frame_nr(idx))
    
    if data_type==1
        image = double(get_pe_new5(dirr, file, frame_nr(idx)));
    elseif data_type==2
        image = permute(double(read_cbf([dirr, file, sprintf('%05d.cbf',frame_nr(idx))]).data),[2 1]);
    else
        disp("Data type not specified properly.")
        break
    end
    
    %     % Background on the border of each image depends both on the absorption
    %     % of the sample (which depends on the shape mostly), and the compton
    %     % background which we assume to be constant
    
    if din_thr > 0
        hotpix=find(image.*(~mask)>din_thr);
        disp([num2str(length(hotpix)),' hot pixel(s)'])
        for hh=(hotpix(:))'
            circ=(xx-xx(hh)).^2+(yy-yy(hh)).^2;
            dinamic_mask1(circ<din_r1^2)=1;      % masking the overexposed pixels + nearby region for 1 image
            dinamic_mask2(circ<din_r2^2)=din_fr; % masking the overexposed pixels for N images
        end
        try
            mask_tot=mask+dinamic_mask1+dinamic_mask2; % updating the mask with the dinamic one
        catch
            mask_tot=mask;
        end
    else
        mask_tot=mask;
    end
    mask_tot=boolean(mask_tot);
    
    image = double(image*myabs(idx)-mybkg);
    % scaling the image by the absorption factor
    % (normilized to the first image) and subtracting
    % the simulated compton scattering

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
    
    % apply P-corection here
    Intens = image ./P;
    %
    final = [pix_h(1:x_size*y_size)', pix_k(1:x_size*y_size)',...
        pix_l(1:x_size*y_size)', Intens(1:x_size*y_size)',...
        mask_tot(1:x_size*y_size)'];
    
    for nn=1:size(final,1)
        if final(nn,1)<n_steps(1) && final(nn,1)>0 && ...
           final(nn,2)<n_steps(2) && final(nn,2)>0 &&  ...
           final(nn,3)<n_steps(3) && final(nn,3)>0 && ...
           final(nn,5)==0
           Int_total(final(nn,1),final(nn,2),final(nn,3)) = ...
               Int_total(final(nn,1),final(nn,2),final(nn,3))+final(nn,4);
           z_total(final(nn,1),final(nn,2),final(nn,3))   = ...
               z_total(final(nn,1),final(nn,2),final(nn,3))+1;
        end
    end
    if din_thr>0
        dinamic_mask1=dinamic_mask1*0; % reinitializing to zero
        idxm=~~dinamic_mask2; % having 1 if masked and 0 if not. Maybe no need for this.
        dinamic_mask2(idxm)=dinamic_mask2(idxm)-1; % decreasing the mask value until it reaches 0 and thus it won't be masked anymore
    end
end

% restoring the TIFF warning
warning('on','MATLAB:imagesci:rtifc:missingPhotometricTag')
end
