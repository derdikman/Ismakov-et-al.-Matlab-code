function [rate_hist,ray_score,ray_angle] = CalculateRayliehAndRateHist(hist_count,hist_time,axis_x)
% Date: 10 of July 2014

condition1 = size(hist_count) == size(hist_time);
condition2 = size(hist_count) == size(axis_x);
condition3 = size(hist_time) == size(axis_x);
if (~((condition1 .* condition2).*condition3))
    disp('ERROR! wrong size of parameters!');
    return;
end

tmp_rate_ang = hist_count./hist_time;
% make NaNs zero!
tmp_rate_ang(find(isnan(tmp_rate_ang))) = 0; %changed from 0 to []
tmp_rate_ang(find(tmp_rate_ang== Inf)) = 0; %same as above

Win=hamming(10);
Win=Win/sum(Win);

%% circular convolution
tmp_rate_ang=cconv(tmp_rate_ang,Win');
rate_hist=tmp_rate_ang((length(Win)/2):length(tmp_rate_ang)-(length(Win)/2));

%find rayleigh score & angle
rate_hist(rate_hist==Inf)=[]; 
norm_val = nansum(rate_hist);
x_vec = nansum(cos(axis_x).*rate_hist);
y_vec = nansum(sin(axis_x).*rate_hist);
vec_len = sqrt(x_vec.^2+y_vec.^2);

ray_score = vec_len/norm_val;
ray_angle = rad2deg(wrapTo2Pi(atan2(y_vec,x_vec)));
