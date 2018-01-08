function [rate_mat]=CreateRateMap(posx,posy,post,spkx,spky,spkt,parms)
%function [rate_mat, spike_mat_smooth]=CreateRateMap(posx,posy,post,spkx,spky,spkt,parms)
max_x = (max(posx)); 
max_y = (max(posy));
min_x = min((posx));
min_y = min((posy));

% divide the environment into spatial bins 
axis_x = min_x-parms.bin_size/2:parms.bin_size:max_x+parms.bin_size/2;
axis_y = min_y-parms.bin_size/2:parms.bin_size:max_y+parms.bin_size/2;
% axis_x = min_x:(max_x-min_x)/50:max_x;
% axis_y = min_y:(max_y-min_y)/50:max_y;
dt=post(2)-post(1);

time_mat = zeros(length(axis_y),length(axis_x));
spike_mat = zeros(length(axis_y),length(axis_x));
rate_mat = zeros(length(axis_y),length(axis_x));

%create time mat (compute how much time the rat spends in each bin)
% find in each moment(time_per_bin) what spatial bin the rat is at and add the time_per_bin to
% 
for i = 1:length(post)
    if ~isnan(posx(i)) && ~isnan(posy(i))
        [min_val,x_ind] =  min(abs(posx(i)-axis_x));
        [min_val,y_ind] =  min(abs(posy(i)-axis_y));
        time_mat(y_ind,x_ind) = time_mat(y_ind,x_ind)+dt;
        
    end
end


%create count mat( count the num of spikes in each bin)
for i = 1:length(spkt)
   if ~isnan(spkx(i)) && ~isnan(spky(i))        %appears to be necessary. remove later if need to-rebekkah
        [min_val,x_ind] =  min(abs(spkx(i)-axis_x));
        [min_val,y_ind] =  min(abs(spky(i)-axis_y));
        spike_mat(y_ind,x_ind)= spike_mat(y_ind,x_ind)+1;
   end
end

%
%spike_mat=Smooth_Rate_Mat(spike_mat,parms);

%time_mat=Smooth_Rate_Mat(time_mat,parms);
% create rate mat

time_mat=SmoothRateMat(time_mat,parms);
spike_mat=SmoothRateMat(spike_mat,parms);

 rate_mat=spike_mat./time_mat;
 rate_mat(rate_mat==inf)=NaN;
 %rate_mat=SmoothRateMat(rate_mat,parms);
 
% spike_mat_smooth =SmoothRateMat(spike_mat, parms);  % want to see spike_mat irrelevant to time 
 
disp('');
 