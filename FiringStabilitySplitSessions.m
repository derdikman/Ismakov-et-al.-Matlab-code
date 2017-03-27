function [all_corrs,all_corrs2, orig_rates_b,orig_rates_e,...
    max_inds_b,max_inds_e,peak_rates_b,peak_rates_e,...
    rate_mats_e,rate_mats_b,zone_mats_begin,zone_mats_end,p,S]= ...
    FiringStabilitySplitSessions(pos_x,pos_y,pos_t,spk_t,...
    zone_mat_total, max_inds,parms)

pos_t_len= length(pos_t);

%finds time indexes of first half and second half of session
if mod(pos_t_len,2) ==0   % even number
    begin_inds= 1:pos_t_len/2;
    end_inds= (pos_t_len/2)+1:pos_t_len;
elseif mod(pos_t_len,2) == 1    % odd number
    begin_inds= 1:(pos_t_len/2)-0.5;
    end_inds= (pos_t_len/2)+1.5:pos_t_len;
end

pos_t_begin= pos_t(begin_inds);   %first half of session
pos_t_end= pos_t(end_inds); %second half of session

spk_t_begin=spk_t(spk_t<pos_t_begin(end)); %spks in first half
spk_t_end=spk_t(spk_t>pos_t_end(1));    %spks in second half

% divide into half
pos_x_begin= pos_x(begin_inds);
pos_y_begin= pos_y(begin_inds);
pos_x_end= pos_x(end_inds);
pos_y_end= pos_y(end_inds);

% build the axis of location when spikes are made for both halves
spk_x_begin=interp1(pos_t_begin,pos_x_begin,spk_t_begin);
spk_y_begin=interp1(pos_t_begin,pos_y_begin,spk_t_begin);
spk_x_end=interp1(pos_t_end,pos_x_end,spk_t_end);
spk_y_end=interp1(pos_t_end,pos_y_end,spk_t_end);

% create rate maps
rate_mat_begin=CreateRateMap(pos_x_begin,pos_y_begin,pos_t_begin,spk_x_begin,spk_y_begin,spk_t_begin,parms);
rate_mat_end=CreateRateMap(pos_x_end,pos_y_end,pos_t_end,spk_x_end,spk_y_end,spk_t_end,parms);

% find gridness score of both halves:


[~,MD_score1,~]=ComputeMovingDirectionalityWithSpeedThreshHold...
    (pos_t_begin,pos_x_begin,pos_y_begin,spk_t_begin,2);

[~,MD_score2,~]=ComputeMovingDirectionalityWithSpeedThreshHold...
    (pos_t_end,pos_x_end,pos_y_end,spk_t_end,2);

if MD_score1 <= 0.15 || MD_score2 <= 0.15
    
    autocorr_1 = Cross_Correlation(rate_mat_begin,rate_mat_begin);
    R_outer_1 = FindROuter(autocorr_1);
    gridness1 = GridnessRadius(autocorr_1,R_outer_1,i);
    autocorr_2 = Cross_Correlation(rate_mat_end,rate_mat_end);
    R_outer_2 = FindROuter(autocorr_2);
    gridness2 = GridnessRadius(autocorr_2,R_outer_2,i);
    
    if gridness1 >= 0.3 || gridness2 >= 0.3
        
        
        rm_size=size(rate_mat_begin);
        nan_sum=sum(sum(rate_mat_begin==0))/(rm_size(1)*rm_size(2));
        if nan_sum < 0.8 % remove  where rat didnt run enough in beginning
            
            % create zone maps
            [zone_mat_b, orig_rates_b, max_inds_b]= getZones(rate_mat_begin);
            [zone_mat_e, orig_rates_e, max_inds_e]= getZones(rate_mat_end);
            
            new_size=size(zone_mat_total);
            old_size=size(zone_mat_b);
            if ~isequal(old_size, new_size)
                [zone_mat_b] = StretchImage(zone_mat_b, size(zone_mat_b), new_size);
            end
            
            old_size=size(zone_mat_e);
            if ~isequal(old_size, new_size)
                [zone_mat_e] = StretchImage(zone_mat_e, size(zone_mat_e), new_size);
            end
            
            % find peak rates of half sessions using full session max inds
            peak_rates_begin= findPeakRates(max_inds,zone_mat_b);
            
            peak_rates_end=findPeakRates(max_inds,zone_mat_e);
            
            peak_rates_b=peak_rates_begin;
            peak_rates_e=peak_rates_end;
            rate_mats_b=rate_mat_begin;
            rate_mats_e=rate_mat_end;
            zone_mats_begin=zone_mat_b;
            zone_mats_end=zone_mat_e;
            
            corr_firing= corrcoef(peak_rates_begin,peak_rates_end);
            
            
            
            
            all_corrs=corr_firing(2);
            
            prb_orig=peak_rates_begin;
            pre_orig=peak_rates_end;
            
            peak_rates_begin(prb_orig==0 | pre_orig==0)= [];
            peak_rates_end(prb_orig==0 | pre_orig==0)= [];
            
            if length(peak_rates_begin) >2
                corr_two= corrcoef(peak_rates_begin,peak_rates_end);
                all_corrs2=corr_two(2);
                results = polyfit(peak_rates_begin,peak_rates_end,1);
                p=results(1);
                S=results(2);
            else
                all_corrs2=nan;
                p=nan;
                S=nan;
            end
            
        else
            all_corrs=nan;
            all_corrs2=nan;
            
            orig_rates_b=nan;
            orig_rates_e=nan;
            max_inds_b=nan;
            max_inds_e=nan;
            peak_rates_b=nan;
            peak_rates_e=nan;
            rate_mats_e=nan;
            rate_mats_b=nan;
            zone_mats_begin=nan;
            zone_mats_end=nan;
            p=nan;
            S=nan;
            
        end
    else
        all_corrs=nan;
        all_corrs2=nan;
        
        orig_rates_b=nan;
        orig_rates_e=nan;
        max_inds_b=nan;
        max_inds_e=nan;
        peak_rates_b=nan;
        peak_rates_e=nan;
        rate_mats_e=nan;
        rate_mats_b=nan;
        zone_mats_begin=nan;
        zone_mats_end=nan;
        p=nan;
        S=nan;
    end
    
else
    all_corrs=nan;
    all_corrs2=nan;
    
    orig_rates_b=nan;
    orig_rates_e=nan;
    max_inds_b=nan;
    max_inds_e=nan;
    peak_rates_b=nan;
    peak_rates_e=nan;
    rate_mats_e=nan;
    rate_mats_b=nan;
    zone_mats_begin=nan;
    zone_mats_end=nan;
    p=nan;
    S=nan;
end

end
