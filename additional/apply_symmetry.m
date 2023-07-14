function [y] = apply_symmetry(x,idx)
%
% [y] = apply_symmetry(x,idx)
% input: unsymmetrised data.
% output: symmetrised data.
% the symmetry operations to be performed are contianed in idx in raws,
% where each axis x,y,z is represented by the numerical code 1,2,3 and 
% flip operations are marked by negative sign.
% for example, the symmetry operations (z,-x,y) and (x,y,-z) correspond
% to idx=[ 3,-1,2; 1,2,-3 ];
%
tic
x(isnan(x))=0;
y=zeros(size(x));
counts=ones(size(x));
for jj = 1:size(idx,1)
    partial = permute(x,abs(idx(jj,:)));
    for ii = 1:size(idx,2)
        if idx(jj,ii) < 0
            partial = flip(partial,ii);
        end
    end
    
    y = y + partial;
    counts = counts + boolean(partial);

    disp(['[',num2str(idx(jj,:)),'] symmetry applied'])
end
y=y./counts;
runtime=toc;
disp([num2str(jj),' symmetry operation applied in ',sprintf('%3.2f',runtime/60),' min'])
end
