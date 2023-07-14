function [data_paf,idx_paf]=paf_karen(data,N,sg,fs,Q,deltaq,delta_shape)
% [data_paf,idx_paf]=paf_karen(data,N,sg,fs,Q,deltaq)
% KAREN based punch anf fill only within a sphere around Braggs;
%
%   N = dimension of the NxNxN moving box; N should be of the order of the
%   characteristic width of the diffuse scattering
%
%   n_steps = diemnsion of the data: n_steps x n_steps x n_steps
%   nominal sg = 3 but larger values might be required for more localized
%   scattering
%   
%   nominal fs = 1 but in this way it might be under filled
%
%   Q = [qh;qk;ql] 3x3 matrix
%   
%   nominal deltaq=0.2 and it's the spherical radius around each Bragg
%   only withing these spheres the punch is possible

% dt=load('Int_total_t_si111_3_flatmono_atten0_dynamicmask.mat');
% data=dt.Int_total./dt.z_total; data(isnan(data))=0; clear dt
% %%
if nargin == 2
    sg=3;
    fs=1;
elseif nargin == 3
    fs=1;
elseif nargin == 5
    deltaq=[1 1 1]*0.2;
    delta_shape=1;
elseif nargin == 6
    delta_shape=1;
end

if length(deltaq)==1
    deltaq=[1 1 1]*deltaq(1);
end

[qk,qh,ql]=meshgrid(Q(2,:),Q(1,:),Q(3,:));

if delta_shape==1
% elipsoid-shape punch
idx_notbraggs=((qh-round(qh)).^2/deltaq(1)^2 + ...
               (qk-round(qk)).^2/deltaq(2)^2 + ...
               (ql-round(ql)).^2/deltaq(3)^2) > 1;
elseif delta_shape==2
    % square-shape punch
    idx1=abs(qh-round(qh))<deltaq(1);
    idx2=abs(qk-round(qk))<deltaq(2);
    idx3=abs(ql-round(ql))<deltaq(3);
    idx_notbraggs=~(idx1.*idx2.*idx3);
end

tic
n_steps=size(data,1);
delta = round(N); % shift of the moving box

if N<=n_steps
    data_paf=data; % initializing;
    data_p=data; % initializing;
    idx_paf=zeros(size(data));

    idx_box=1:N;
    for ii=1:delta:n_steps-N+1
        for jj=1:delta:n_steps-N+1
            for kk=1:delta:n_steps-N+1
%                 [aa,bb,cc]=meshgrid(idx_box+(ii-1),idx_box+(jj-1),idx_box+(kk-1));
                mbox=data(idx_box+(ii-1),idx_box+(jj-1),idx_box+(kk-1)); idxf = ~~mbox; % zero values should not be considered in the proccess
                % *****************************************************
                % idx_box+(ii-1),idx_box+(jj-1),idx_box+(kk-1) are the
                % indices of the current box within the whole dataset
                % *****************************************************
                if find(idxf)>0
                    medt = median(mbox(idxf)); % median of the data in the moving box
%                         disp(['median = ',num2str(medt)])
                    mad = median(abs(mbox(idxf) - medt)); % differential median for the threshold determination
                    thr = medt + sg * 1.4826 * mad; % threshold for the outliers
%                     thr = medt + 3 * 1.4826 * mad; % threshold for the outliers
%                         disp(['threshold = ',num2str(thr)])
                    fill = medt + 2.2 * mad * fs; % filling value for the punched pixels
%                     fill = medt + 2.2 * mad*5; % filling value for the punched pixels
%                         disp(['filling = ',num2str(fill)])
                    if mad>0
                        idx_out = find(mbox > thr);
                        % ************************************************
                        % idx_out are the linear indices of the pixels to
                        % be punched
                        % ************************************************
                    else
                        idx_out=0;
                    end
                else
                    idx_out=0;
                end
                
                if idx_out ~= 0
                    if ~isempty(idx_out)
                        % updating idx_out with only Bragg positions
                        check_nobraggs=idx_notbraggs(idx_box+(ii-1),idx_box+(jj-1),idx_box+(kk-1));
                        i_nobragg=find(check_nobraggs(idx_out));
                        %                 i_nobragg=[];
                        %                 for i_check=1:length(idx_out)
                        %                     if find(idx_notbraggs == idx_out(i_check))
                        %                         i_nobragg=[i_nobragg;idx_out(i_check)];
                        %                     end
                        %                 end
                        if ~isempty(i_nobragg)
                            idx_out(i_nobragg)=[];
                        end
                        if ~isempty(idx_out)
                            % punch and fill
                            mbox(idx_out)=0;
                            data_p(idx_box+(ii-1),idx_box+(jj-1),idx_box+(kk-1))=mbox; % only punched data
                            mbox(idx_out)=fill;
                            data_paf(idx_box+(ii-1),idx_box+(jj-1),idx_box+(kk-1))=mbox; % punched and filled data (KAREN)
                            mbox(idx_out)=nan;
                            idx_paf(idx_box+(ii-1),idx_box+(jj-1),idx_box+(kk-1))=mbox; % indices ofthe punched pixels
                        end
                    end
                end
            end
        end
%         disp([ii jj kk])
    end

    paftime=toc;
    disp(['punch & fill (KAREN) in ',num2str(paftime/60),' min'])
else
    disp('ERROR - the moving-box size is larger than the data')
    disp(['data are: [',num2str(size(data)),'] px^3'])
end
idx_paf=find(isnan(idx_paf)); % linear indices of the punched data
%%
end
