cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis');
load('3sets G3MD15PF3 data and results');

%info for rate map creation
parms.bin_size= 3;
parms.sigma= 1.5;

len= length(rate_mats_all); % length of # of cells

%initialize peak rates of beginning (b) and end (e) of session
peak_rates_b=cell(1,len);
peak_rates_e=cell(1,len);

twoD_corrs=nan(1,len);
twoD_zone_corrs=nan(1,len);
all_corrs=nan(1,len);
all_corrs2=nan(1,len);

rate_mats_b=cell(1,len);
rate_mats_e=cell(1,len);
zone_mats_begin=cell(1,len);
zone_mats_end=cell(1,len);
orig_rates_b=cell(1,len);
orig_rates_e=cell(1,len);
max_inds_b=cell(1,len);
max_inds_e=cell(1,len);

for i =1:len
    
    if i ~= 69 
        
    pos_x=pos_x_all{i};
    pos_y=pos_y_all{i};
    pos_t=pos_t_all{i};
    spk_x=spk_x_all{i};
    spk_y=spk_y_all{i};
    spk_t=spk_t_all{i};
    
    rate_mat_total=rate_mats_all{i};
    zone_mat_total=zone_mats_all{i};
    max_inds=max_indices{i};
  
   [all_corrs(i),all_corrs2(i), orig_rates_b{i},orig_rates_e{i},...
    max_inds_b{i},max_inds_e{i},peak_rates_b{i},peak_rates_e{i},...
    rate_mats_e{i},rate_mats_b{i},zone_mats_begin{i},zone_mats_end{i},all_p(i),all_s(i)]= ...
    FiringStabilitySplitSessions(pos_x,pos_y,pos_t,spk_t,...
    zone_mat_total, max_inds,parms);
     
    i
    end
end % of # of cells

save('half sessions G3MD15PF3 info results diff criteria', 'peak_rates_b', 'peak_rates_e',...
    'rate_mats_e', 'rate_mats_b', 'zone_mats_begin', 'zone_mats_end',...
    'all_corrs', 'all_corrs2','twoD_zone_corrs', ...
    'max_inds_b','max_inds_e','orig_rates_e','orig_rates_b','all_p','all_s');

