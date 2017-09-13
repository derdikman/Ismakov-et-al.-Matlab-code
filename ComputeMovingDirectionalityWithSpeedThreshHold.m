function [rate_hist,ray_score,ray_angle]=ComputeMovingDirectionalityWithSpeedThreshHold...
    (pos_t,pos_x,pos_y,spk_t,speed_thresh_hold)


dt=median(diff(pos_t));
%speed = sqrt(diff(pos_x).^2 + diff(pos_y).^2) / dt;

% pos_x=pos_x/5;
% pos_y=pos_y/5;

speed = sqrt(diff(pos_x).^2 + diff(pos_y).^2) / dt;


% TODO smothing - or to do the smothing outside of the function (on pos_x &
% pos_y)

% the moving direction of the animal through out the trail
win=hamming(11);
win=win/sum(win);

diff_x=conv(diff(pos_x),win,'same');
diff_y=conv(diff(pos_y),win,'same');

time_phi= wrapTo2Pi(atan2(diff_y,diff_x));

%time_phi_vec(1:length(time_phi),1)=time_phi;
spk_phi=interp1(pos_t,[0 time_phi'],spk_t);
spk_speed = interp1(pos_t,[0 speed'],spk_t);
% the head direction of the animal when spike has ocurred
%count_phi = atan2(diff(spk_x), diff(spk_y));
% cut all the data under the thresh hold
time_phi = time_phi(find(speed>=speed_thresh_hold));
count_phi = spk_phi(find(spk_speed>=speed_thresh_hold));

step_len=deg2rad(3);
ang_ax =0:step_len:2*pi;

count = hist(count_phi,ang_ax);
time = hist(time_phi,ang_ax)*dt;

[rate_hist,ray_score,ray_angle] = CalculateRayliehAndRateHist(count,time,ang_ax);

