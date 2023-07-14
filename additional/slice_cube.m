function [datacut,xa,ya,za]=slice_cube(data,vx,vy,vz,dx,dy,dz,bin)

switch nargin
    case 5
        dy=dx; dz=dx;
        bin=0;
    case 7
        bin=0;
end

if numel(dx) == 2
    i1=find(vx<=dx(1),1,'last'); i2=find(vx>=dx(2),1,'first');
    if isempty(i1); i1 = 1; end; if isempty(i2); i2 = length(vx); end
    idxx=i1:i2;
else
    idxx=1:numel(vx);
end
if numel(dy) == 2
    i1=find(vy<=dy(1),1,'last'); i2=find(vy>=dy(2),1,'first');
    if isempty(i1); i1 = 1; end; if isempty(i2); i2 = length(vy); end
    idxy=i1:i2;
else
    idxy=1:numel(vy);
end
if numel(dz) == 2
    i1=find(vz<=dz(1),1,'last'); i2=find(vz>=dz(2),1,'first');
    if isempty(i1); i1 = 1; end; if isempty(i2); i2 = length(vz); end
    idxz=i1:i2;
else
    idxz=1:numel(vz);
end

xa=vx(idxx); ya=vy(idxy); za=vz(idxz);
datacut = data(idxx,idxy,idxz);  % intensities
% dt      = zdata(idxx,idxy,idxz); % number of counts
dt = datacut ~= 0;
db=find(bin==1);
if numel(bin) == 3
    switch db
        case 1
            datacut=squeeze(sum(datacut,1))./squeeze(sum(dt,1));
            xa=ya; ya=za; za=mean(dx);
        case 2
            datacut=squeeze(sum(datacut,2))./squeeze(sum(dt,2));
            ya=za; za=mean(dy);
        case 3
            datacut=squeeze(sum(datacut,3))./squeeze(sum(dt,3));
            za=mean(dz);
    end
end
datacut(isnan(datacut))=0;
end