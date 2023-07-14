function [pos1,wid1] = gauss_fit_b(x,y)
    
% starting parameters (do a numerical analysis):
% amp1: maximum - background (see below)
[amp1_start,index_start]=max(y);

% pos1: x(maximum)
pos1_start=x(index_start);

% bg: (av(left3points)+av(right3points))/2
left3points=(y(1)+y(2)+y(3))/3;
right3points=(y(length(y))+y(length(y)-1)+y(length(y)-2))/3;
bg_start=double((left3points+right3points)/2);

% slope: (av(left3points)-av(right3points))/x-range    
xrange=(x(length(x))-x(1));
slope_start=double((right3points-left3points)/xrange);

% amp1: maximum - background
amp1_start=amp1_start-bg_start;

nobg_y=y-bg_start-slope_start*(x-pos1_start);

% right half point
p=0;
while nobg_y(index_start+p) > amp1_start/2
    p=p+1;
end
y_dist=-nobg_y(index_start+p)+nobg_y(index_start+p-1);
half_dist=-nobg_y(index_start+p)+amp1_start/2;
factor=half_dist/y_dist;
x_dist=x(index_start+p)-x(index_start+p-1);
rhp=x(index_start+p)-factor*x_dist;

% left half point
q=0;
while nobg_y(index_start-q) > amp1_start/2
    q=q+1;
end
y_dist=nobg_y(index_start-q)-nobg_y(index_start-q+1);
half_dist=nobg_y(index_start-q)-amp1_start/2;
factor=half_dist/y_dist;
x_dist=-x(index_start-q)+x(index_start-q+1);
lhp=x(index_start-q)+factor*x_dist;

fwhm_start=rhp-lhp;
midpoint=(rhp+lhp)/2;
wid1_start=fwhm_start/(2*0.693);
pos1_start=midpoint;

amp1=amp1_start;
pos1=pos1_start;
wid1=wid1_start;
bg=bg_start;
slope=slope_start;

integral=amp1*wid1;

end
