function [det_dist,orgx,orgy,hor_angle,vert_angle]=poni2org(fullpath,dx)
% [det_dist,orgx,orgy,hor_angle,vert_angle]=poni2org(fullpath,dx)
% fullpath is with the filename and extension
% dx = is the pixel size in mm
% O. Ivashko 2021

% identifiers
dist_id  = 'Distance: ';
poni1_id = 'Poni1: ';
poni2_id = 'Poni2: ';
rot1_id  = 'Rot1: ';
rot2_id  = 'Rot2: ';

% initialization
num=size(fullpath,1); % number of distinct files
d=0;
poni1=0;
poni2=0;
rot1=0;
rot2=0;
det_dist=zeros(1,num);
orgx=zeros(1,num);
orgy=zeros(1,num);
hor_angle=zeros(1,num);
vert_angle=zeros(1,num);

for idx=1:num
    % loading ponifile and extracting values
    fid = fopen(fullpath(idx,:));
    if fid==-1
        disp('File not found')
        disp(fullpath(idx,:))
    else
    dataline=fgetl(fid);
    ok=1;
    while ok
        if strncmp(dataline,dist_id,length(dist_id))
            ok=0;
            d = str2double(dataline(length(dist_id)+1:end));
        else
            dataline=fgetl(fid);
        end
    end

    dataline=fgetl(fid);
    if strncmp(dataline,poni1_id,length(poni1_id))
        poni1 = str2double(dataline(length(poni1_id)+1:end));
    end

    dataline=fgetl(fid);
    if strncmp(dataline,poni2_id,length(poni2_id))
        poni2 = str2double(dataline(length(poni2_id)+1:end));
    end

    dataline=fgetl(fid);
    if strncmp(dataline,rot1_id,length(rot1_id))
        rot1 = str2double(dataline(length(rot1_id)+1:end));
    end

    dataline=fgetl(fid);
    if strncmp(dataline,rot2_id,length(rot2_id))
        rot2 = str2double(dataline(length(rot2_id)+1:end));
    end
    end
    fclose(fid);

    % converting to ORG reference
    det_dist(idx) = d*1e3 / cos(rot1) / cos(rot2);

    % in pyFAI images are flipped along the horizontal axis
    % (0,0) point is still the same as in Matlab
    % this results in change of sign of rot1 (rotation around the vertical)
    % rot2 (rotation around the horizontal) instead is unchanged
    orgx(idx) = (poni1*1e3 + d*1e3 * tan( rot2)) / dx;
    orgy(idx) = (poni2*1e3 + d*1e3 * tan(-rot1)) / dx;

    hor_angle(idx)  = rot2;
    vert_angle(idx) =-rot1;
end

if num==1
    disp('###################################')
    disp('### conversion from PONI to ORG ###')
    disp('###################################')
    disp([sprintf('%s',pad('det_dist',10,'right')),' = ',sprintf('%3.4f',det_dist),';     % mm'])
    disp([sprintf('%s',pad('orgx',10,'right')),' = ',sprintf('%3.4f',orgx),';     % px'])
    disp([sprintf('%s',pad('orgy',10,'right')),' = ',sprintf('%3.4f',orgy),';     % px'])
    disp([sprintf('%s',pad('hor_angle',10,'right')),' = ',sprintf('%3.8f',hor_angle),';  % rad'])
    disp([sprintf('%s',pad('vert_angle',10,'right')),' = ',sprintf('%3.8f',vert_angle),';  % rad'])
    disp('###################################')
    disp('###################################')
else
        disp('###################################')
    disp('### conversion from PONI to ORG ###')
    disp('###################################')
    disp([sprintf('%s',pad('det_dist',10,'right')),' = [',sprintf('%3.4f ',det_dist),'];   % mm'])
    disp([sprintf('%s',pad('orgx',10,'right')),' = [',sprintf('%3.4f ',orgx),'];   % px'])
    disp([sprintf('%s',pad('orgy',10,'right')),' = [',sprintf('%3.4f ',orgy),'];   % px'])
    disp([sprintf('%s',pad('hor_angle',10,'right')),' = [',sprintf('%3.8f ',hor_angle),'];   % rad'])
    disp([sprintf('%s',pad('vert_angle',10,'right')),' = [',sprintf('%3.8f ',vert_angle),'];   % rad'])
    disp('###################################')
    disp('###################################')
end

end