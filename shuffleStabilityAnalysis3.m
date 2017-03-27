function shuffleStabilityAnalysis3

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis')
load('half sessions G3MD15PF3 info results diff criteria.mat',...
    'orig_rates_b','orig_rates_e','max_inds_b','max_inds_e', 'all_p','all_s')

load('3sets G3MD15PF3 data and results.mat',...
    'zone_mats_all','max_indices', 'PF_radii')

shuffle_times=500;

mean_corrs=nan(1,shuffle_times);
mean_corrs2=nan(1,shuffle_times);
mean_2D_corrs=nan(1,shuffle_times);

for k= 1:shuffle_times
    
    len=length(zone_mats_all);
    
    all_corrs=nan(1,len);
    all_corrs2=nan(1,len);
    twoD_zone_corrs=nan(1,len);
    
    for i=1:len
        
        PF_radius=PF_radii(i);
        zone_mat_total=zone_mats_all{i};
        max_inds=max_indices{i};
        
        max_inds_end=max_inds_e{i};
        max_inds_begin=max_inds_b{i};
        rates_b=orig_rates_b{i};
        rates_e=orig_rates_e{i};
        
        if ~isnan(max_inds_begin)
        
        %shuffle rates
        rates_b=Shuffle(rates_b);
        rates_e=Shuffle(rates_e);
        
        size_zm=size(zone_mat_total);
        
        [zone_mat_e,~]=CreateZoneMat(size_zm, PF_radius, max_inds_end, rates_e);
        [zone_mat_b,~]=CreateZoneMat(size_zm, PF_radius, max_inds_begin, rates_b);
        
        % find peak rates of half sessions using full session max inds
       peak_rates_begin=findPeakRates(max_inds,zone_mat_b);
       peak_rates_end=findPeakRates(max_inds,zone_mat_e); 
        
        % find different correlations
        corr_firing= corrcoef(peak_rates_begin,peak_rates_end);
        all_corrs(i)=corr_firing(2);
        
        prb_orig=peak_rates_begin;
        pre_orig=peak_rates_end;
        
        peak_rates_begin(prb_orig==0 | pre_orig==0)= [];
        peak_rates_end(prb_orig==0 | pre_orig==0)= [];
        
        if length(peak_rates_begin) >2
            corr_two= corrcoef(peak_rates_begin,peak_rates_end);
            all_corrs2(i)=corr_two(2);
            results = polyfit(peak_rates_begin,peak_rates_end,1);
            all_sh_p(i)=results(1);
            all_sh_s(i)= results(2);
        else
            all_corrs2(i)=nan;
            all_sh_p(i)=nan;
            all_sh_s(i)=nan;
        end
        
        twoD_zone_corrs(i)=corr2(zone_mat_b,zone_mat_e);
        
        end 
    end
    
    mean_corrs(k)=nanmean(all_corrs);
    mean_corrs2(k)=nanmean(all_corrs2);
    mean_2D_corrs(k)=nanmean(twoD_zone_corrs);
    mean_p(k) = nanmean(all_sh_p);
    mean_s(k) = nanmean(all_sh_s);
    
    k
    
    save('shuffled split session firing stability analysis diff criteria', ...
    'mean_corrs','mean_corrs2','mean_2D_corrs','mean_p','mean_s');

end


% actual mean corr= 0.4064
% actual mean corr2= 0.4184 0.4938
% actual 2D corr= 0.4175

figure; hist(mean_p)
figure;hist(mean_s)
