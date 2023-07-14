function [img] = get_pe_new5(dirr, file, filenumber)
% [img] = get_pe_new5(dirr, file, filenumber)

% disabling the TIFF warning to increase the speed
warning('off','MATLAB:imagesci:rtifc:missingPhotometricTag')

%disp(filenumber);
if nargin==0
    disp('Filenumber not given')
    return
else
    filename = [dirr,file,sprintf('%05d',filenumber),'.tif'];
end
% disp(filename)
img=imread(filename);
% img=imread_noerror(filename);

% % enabling the TIFF warning to increase the speed
% warning('on','MATLAB:imagesci:rtifc:missingPhotometricTag')
end