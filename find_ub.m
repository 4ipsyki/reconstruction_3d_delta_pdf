base_dirr        = '/asap3/petra3/gpfs/p07/2019/data/11006746/';
name_ip          = 'ip_700_allpos_02';
pars_dirr        = [base_dirr,'processed/inputpar/'];

addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/core/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/additional/']);
addpath([base_dirr,'processed/reconstructions/ybco_oleh/scripts/find_ub_matrix/']);
load([pars_dirr,name_ip,'.mat']);

nr=5;
    dirrname=dirr{nr}; flname=file{nr};
    ox=orgx(nr);oy=orgy(nr);
    sdd=det_dist(nr);
    h_angle=hor_angle(nr);
    v_angle=vert_angle(nr);
    om_pos=rot_angle{nr};
    frames=frame_nr{nr};
    clear mask; load([save_dirr,'masks_700.mat']);
    mask=mask{nr};

%% loading peak-angle matrix or loading frames and finding peaks
try
    load([save_dirr,flname,'list_xpix_ypix_angle.mat']);
catch

    list_xpix_ypix_angle=[0,0,0];
    nsum=100; % number of images to be summed to search for a peak, i.e. moving average size
    thres=10.0; % threshold for peak detection
    niter=round(n_pics/nsum); % number of loops
    
    for nb=1:niter
        image_ary=zeros(x_size,y_size,nsum);
        filenumber=frames((1:nsum)+nsum*(nb-1));
        parfor n=1:length(filenumber)
            display(filenumber(n));
            img = single(get_pe_new5(dirrname, flname, filenumber(n)));
            image_ary(:,:,n) = img;
        end

        show_img = sum(image_ary,3)/nsum;
        show_img(boolean(mask))=0;

        % find peaks on the summed image and display them
        pee=FastPeakFind(show_img, thres);
        imagesc(show_img), shading flat, colormap viridis
        caxis([0 thres+2])
        hold on
        plot(pee(1:2:end),pee(2:2:end),'r+')
        hold off
        pause(0.2);

        % determine x_pixel, y_pixel, rot_angle for all peaks
        for o=1:2:size(pee,1)
            m=pee(o);
            n=pee(o+1);
            scan=squeeze(image_ary(n,m,:));
            if max(scan) > 0
                if max(scan) > scan(1)*2+scan(size(scan,1))*2
                    try
                        [pos,fwhm] = gauss_fit_b(om_pos((1:nsum)+nsum*(nb-1)),scan);
                    catch
                        disp("couldn't fit a peak")
                    end
                    if fwhm<3
                        list_xpix_ypix_angle=[list_xpix_ypix_angle;n,m,pos];
                    end
                end
            end
        end
    end
    save([save_dirr,flname,'list_xpix_ypix_angle.mat'],'list_xpix_ypix_angle');
end
%% transform in rez. space with unit matrix
pars = [ox,oy,dx,dy,sdd,lambda,ki_vec,det_x,det_y,...
        h_angle,v_angle];
abc=[1,0,0;  0,1,0;  0,0,1];
list_as_bs_cs=zeros(size(list_xpix_ypix_angle));
for peak=2:size(list_xpix_ypix_angle,1)
    [a_star,b_star,c_star] = ...
        hkl_from_abc(list_xpix_ypix_angle(peak,1),...
                                  list_xpix_ypix_angle(peak,2),...
                                  list_xpix_ypix_angle(peak,3),...
                                  abc,pars);
    list_as_bs_cs(peak+1,:)=[a_star,b_star,c_star];
end

% plotting list_as_bs_cs. when rotating the pekas should form a regular pattern
% if this is not the case, some of the input parameters are not correct
list_as_bs_cs=list_as_bs_cs(2:size(list_as_bs_cs,1),:,:,:);
scatter3(list_as_bs_cs(:,1),list_as_bs_cs(:,2),list_as_bs_cs(:,3))
%view([-183.5 7.18])
xlabel('x');
ylabel('y');
zlabel('z');

%% powder-like plot to select the desired Q-value for ub-matrix determination
% determine the distance between all peaks
ll=length(list_as_bs_cs);
myn=(ll-1)/2; % number of rotations to do; if mod is 0.5 it means that the
              % last subtraction vector is twice reduntand.
              % in this case only the first half has to be considered
if round(myn)-myn == 0.5
    maxshift = myn-0.5;
    additionalshift = 1;
else
    maxshift = myn;
    additionalshift = 0;
end
list_length=zeros(ll*myn,4);
for shift=1:maxshift
    diff_coord=list_as_bs_cs-circshift(list_as_bs_cs,shift);
    dist=vecnorm(diff_coord,2,2);
    list_length((1:ll)+ll*(shift-1),:)=[dist,diff_coord];
end
if additionalshift
    shift=shift+1;
    diff_coord=list_as_bs_cs-circshift(list_as_bs_cs,shift);
    dist=vecnorm(diff_coord,2,2);
    list_length((1:ll/2)+ll*(shift-1),:)=[dist(1:ll/2),diff_coord(1:ll/2,:)];
end

% plot histogram of distances
% show chosen lattice parameters in red, green and blue
myhist=histogram(list_length(:,1),0:0.001:1);
hh=1; % a*-axis peak in r.l.u.
kk=1; % b*-axis peak in r.l.u.
ll=1; % c*-axis peak in r.l.u.
% it's better to choose the lowest-deltaQ peaks possible.
% note that even if the (1,0,0) peak is not allowed, in the differential map
% this peak is usually present.\

% range for the peak selection
delta=0.06;
hold on
mymax=max(myhist.Values);
plot( [(hh/lattice(1))*(1-delta),(hh/(lattice(1)))*(1-delta)],[0 mymax],'red')
plot( [(hh/lattice(1))*(1+delta),(hh/(lattice(1)))*(1+delta)],[0 mymax],'red')
plot( [(kk/lattice(2))*(1-delta),(kk/lattice(2))*(1-delta)],[0 mymax],'green')
plot( [(kk/lattice(2))*(1+delta),(kk/lattice(2))*(1+delta)],[0 mymax],'green')
plot( [(ll/lattice(3))*(1-delta),(ll/lattice(3))*(1-delta)],[0 mymax],'blue')
plot( [(ll/lattice(3))*(1+delta),(ll/lattice(3))*(1+delta)],[0 mymax],'blue')
hold off
xlim([0 1])

%% plot of chosen-Q-range peaks
% show position of chosen peaks in deltaQ-space
test_a=list_length((list_length(:,1)>(hh/lattice(1))*(1-delta) & list_length(:,1)<(hh/lattice(1))*(1+delta)),:,:,:);
test_b=list_length((list_length(:,1)>(kk/lattice(2))*(1-delta) & list_length(:,1)<(kk/lattice(2))*(1+delta)),:,:,:);
test_c=list_length((list_length(:,1)>(ll/lattice(3))*(1-delta) & list_length(:,1)<(ll/lattice(3))*(1+delta)),:,:,:);


clf, hold on
scatter3(test_a(:,2),test_a(:,3),test_a(:,4),'r','sizedata',100);
scatter3(test_b(:,2),test_b(:,3),test_b(:,4),'g','sizedata',100)
scatter3(test_c(:,2),test_c(:,3),test_c(:,4),'b','sizedata',100);
legend('q_h','q_k','q_l')

% the approximate amount of clusters in the figure needs to be choosen manually
n_cluster=20;
    cluster_a = clusterdata([test_a(:,2),test_a(:,3),test_a(:,4)],n_cluster);
    figure(777)
    scatter3(test_a(:,2),test_a(:,3),test_a(:,4),10,cluster_a); hold on

    cluster_b = clusterdata([test_b(:,2),test_b(:,3),test_b(:,4)],n_cluster);
    figure(777)
    scatter3(test_b(:,2),test_b(:,3),test_b(:,4),10,cluster_b)

    cluster_c = clusterdata([test_c(:,2),test_c(:,3),test_c(:,4)],n_cluster);
    figure(777)
    scatter3(test_c(:,2),test_c(:,3),test_c(:,4),10,cluster_c); hold off 

%% making clusters. calcutaling angles between clusters. selecting the closes to nominal with larger population. determining ub-matrix
figure(22)
clf, dime(10,20)
subplot(3,1,1)
my2hist_a=histogram(cluster_a);
title('cluster a')

% choosing the threshold for the selection of clusters to be considered for the ub matrix
cluster_size_min=100;

% position of all clusters with more than cluster_size_min
huhu=[test_a,cluster_a];
pos_clusters_a=[0,0,0];
for i=1:length(my2hist_a.Values)
    if my2hist_a.Values(i)>cluster_size_min
        haufen2=huhu(huhu(:,5)==i,:,:,:,:,:);
        pos_clusters_a=[pos_clusters_a;mean(haufen2(:,2)),mean(haufen2(:,3)),mean(haufen2(:,4))];
    end
end

figure(22)
subplot(3,1,2)
my2hist_b=histogram(cluster_b);
title('cluster b')

% position of all clusters with more than cluster_size_min
huhu=[test_b,cluster_b];
pos_clusters_b=[0,0,0];
for i=1:length(my2hist_b.Values)
    if my2hist_b.Values(i)>cluster_size_min
        haufen2=huhu(huhu(:,5)==i,:,:,:,:,:);
        pos_clusters_b=[pos_clusters_b;mean(haufen2(:,2)),mean(haufen2(:,3)),mean(haufen2(:,4))];
    end
end

figure(22)
% clf
subplot(3,1,3)
my2hist_c=histogram(cluster_c);
title('cluster c')

% position of all clusters with more than cluster_size_min
huhu=[test_c,cluster_c];
pos_clusters_c=[0,0,0];
for i=1:length(my2hist_c.Values)
    if my2hist_c.Values(i)>cluster_size_min
        haufen2=huhu(huhu(:,5)==i,:,:,:,:,:);
        pos_clusters_c=[pos_clusters_c;mean(haufen2(:,2)),mean(haufen2(:,3)),mean(haufen2(:,4))];
    end
end

% find 3-tupels of clusters whose angles fits best to chosen angle
% angles between all clusters
angles_abc=[0,0,0, 0,0,0, 0,0,0, 0,0,0];
% c-a angles
for j=2:size(pos_clusters_c,1)
    for k=2:size(pos_clusters_b,1)
    %for j=2
        for l=j:size(pos_clusters_a,1)
        % angle between the haufen
        u=pos_clusters_c(j,:);
        v=pos_clusters_b(k,:);
        w=pos_clusters_a(l,:);
        angle_wv = atan2d(norm(cross(w,v)),dot(w,v));
        angle_uw = atan2d(norm(cross(u,w)),dot(u,w));
        angle_uv = atan2d(norm(cross(u,v)),dot(u,v));
        angles_abc=[angles_abc;u,v,w,angle_wv,angle_uw,angle_uv];
        end
    end
end
huhu=angles_abc
% select best tupel:
best=(angles_abc(:,10)-lattice(4)).^2 + (angles_abc(:,11)-lattice(5)).^2 + (angles_abc(:,12)-lattice(6)).^2;
[M,I]=min(best);

ax1 = [angles_abc(I,1,:), angles_abc(I,2,:), angles_abc(I,3,:)]/ll;   % a
ax2 = [angles_abc(I,4,:), angles_abc(I,5,:), angles_abc(I,6,:)]/kk;
ax3 = [angles_abc(I,7,:), angles_abc(I,8,:), angles_abc(I,9,:)]/hh;

len = [norm(ax1), norm(ax2), norm(ax3)]
win = [atan2d(norm(cross(ax1,ax2)),dot(ax1,ax2)), atan2d(norm(cross(ax2,ax3)),dot(ax2,ax3)), atan2d(norm(cross(ax1,ax3)),dot(ax1,ax3))]

% initial matrix 
vol=dot(ax1,cross(ax2,ax3));
a_axis=cross(ax2,ax3)/vol;
b_axis=cross(ax3,ax1)/vol;
c_axis=cross(ax1,ax2)/vol;
abc=[a_axis;b_axis;c_axis]'; % tentative ub matrix

% % ordering axes in size and reshaping the ub matrix accordingly
% nra=norm(a_axis); nrb=norm(b_axis); nrc=norm(c_axis);
% [axs_length,order]=sort([nra,nrb,nrc])
% abc=abc(:,order)

n_cluster=50;
cluster_a = clusterdata([test_a(:,2),test_a(:,3),test_a(:,4)],n_cluster);
cluster_b = clusterdata([test_b(:,2),test_b(:,3),test_b(:,4)],n_cluster);
cluster_c = clusterdata([test_c(:,2),test_c(:,3),test_c(:,4)],n_cluster);

figure(23)
hold on
if exist('pl','var');delete(pl);end
pl(1)=plot3([0,ax1(1)]*ll,[0,ax1(2)]*ll,[0,ax1(3)]*ll,'r','linewidth',2);
pl(2)=plot3([0,ax2(1)]*kk,[0,ax2(2)]*kk,[0,ax2(3)]*kk,'g','linewidth',2);
pl(3)=plot3([0,ax3(1)]*hh,[0,ax3(2)]*hh,[0,ax3(3)]*hh,'b','linewidth',2);
hold off, axis equal
xlabel('x'),ylabel('y'),zlabel('z')

%% indexing the found peaks
delete newlist.dat;
list_hkl=[0,0,0];

myrvalue=0;
for peak=1:size(list_xpix_ypix_angle,1)
    [a_star,b_star,c_star]=hkl_from_abc(list_xpix_ypix_angle(peak,1),...
                                           list_xpix_ypix_angle(peak,2),...
                                           list_xpix_ypix_angle(peak,3),...
                                           abc, pars);
    list_hkl=[list_hkl;a_star,b_star,c_star];
    myrvalue=myrvalue+(abs(a_star-round(a_star))+...
                       abs(b_star-round(b_star))+...
                       abs(c_star-round(c_star)))/3;
    % newlist: m,n,image_nr, h,k,l
    reflection=[list_xpix_ypix_angle(peak,1),list_xpix_ypix_angle(peak,2),...
        list_xpix_ypix_angle(peak,3),round(a_star),round(b_star),round(c_star)];
    % write within the selected range (around integer values) to a list
    myrange = 0.1;
    if peak > 1 && abs((a_star-round(a_star))/a_star) < myrange &&...
                   abs((b_star-round(b_star))/b_star) < myrange &&...
                   abs((c_star-round(c_star))/c_star) < myrange
        fileID = fopen('newlist.dat','a');
        fmt = '%5d %5d %5d %5d %5d %5d\n';
        fprintf(fileID,fmt,reflection);
        fclose(fileID);
    end
end
hihi=list_hkl
disp('residual:')
myrvalue=myrvalue/size(list_xpix_ypix_angle,1)
%% refine the initial matrix abc from abov by using all reflections of the list list_xpix_ypix_imnum
x0= double(abc);

% A*x ≤ b
A=[];
b=[];
% Aeq*x = beq
Aeq=[];
beq=[];
lb=[];
ub=[];

% saving parameters that are used in matrix_constraint.m
pars = [ox,oy,dx,dy,sdd,lambda,ki_vec,det_x,det_y,h_angle,v_angle,lattice];
save('pars.mat','pars')

% optimizing the ub matrix
nonlcon = @matrix_constraint;
x = fmincon(@matrix_fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
% x is the new optimized ub matrix
% define abc=x and repeat the previous section about the peak indexing:
% 1) see if the r-value improved; 2) repeat the optimization
% repeat until x is not changing.
% update the "ub" variable in the input-parameter file.

% %% OPTIONAL: refining detector orientation
% % when the peak indexing is not working properly, one can try to optimize
% % the detector parameters.
% % NOTE that it's linked to the current ub matrix ("x" variable)!
% % This means that the tentative ub matrix has to be almost correct.
% abc0=x;
% 
% % A*x ≤ b
% A=[];
% b=[];
% % Aeq*x = beq
% Aeq=[];
% beq=[];
% lb=[];
% ub0=[];
% 
% % saving parameters used in matrix_constraint_detector_orientation.m
% pars_do = [dx,dy,lambda,ki_vec,det_x,det_y,lattice,[abc0(1,:),abc0(2,:),abc0(3,:)]];
% save('pars_detector_orientation.mat','pars_do');
% 
% xdo0=[ox,oy,sdd,h_angle,v_angle];
% 
% addpath([base_dirr,'processed/scripts/'])
% nonlcon = @matrix_constraint_detector_orientation;
% xdo = fmincon(@matrix_fun_detector_orientation,xdo0,A,b,Aeq,beq,lb,ub0,nonlcon);
% 
% ox=xdo(1); oy=xdo(2); sdd=xdo(3); h_angle=xdo(4); v_angle=xdo(5);
% pars = [ox,oy,dx,dy,sdd,lambda,ki_vec,det_x,det_y,h_angle,v_angle,lattice];
% % iterate with new "x" refinement, and peak-list section
% % if the new parameters are reasobale, update the input-parameter file accordingly.


