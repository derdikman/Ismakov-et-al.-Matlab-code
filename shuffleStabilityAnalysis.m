function shuffleStabilityAnalysis

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis')
% load('same arenas corr coef results.mat', ...
%     'all_zm_2','PF_rad_all',...
%     'orig_rates_1','orig_rates_2', 'all_max_inds_1','all_max_inds_2')
%  load('corr coef results of rescaled arenas.mat', ...
 load('corr coef results of same arenas COMBINED.mat', ...
     'all_zm_2',...
     'orig_rates_1','orig_rates_2', 'all_max_inds_1','all_max_inds_2',...
     'PF_radii2')

load('shuffled same areans firing stability analysis.mat')
shuffle_times=10000;

%mean_corr=nan(1,shuffle_times);
%mean_corr2=nan(1,shuffle_times);
for k=8210:shuffle_times
    
    len=length(all_zm_2);
    
    all_corrs=nan(1,len);
    all_corrs2=nan(1,len);
    for i=1:len
        
        zm2=all_zm_2{i};
        rates_1=orig_rates_1{i};
        orig_rates2=orig_rates_2{i};
        max_inds1=all_max_inds_1{i};
        max_inds2=all_max_inds_2{i};
        PF_rad=PF_radii2(i);
       
        orig_rates2=Shuffle(orig_rates2);
        
       
        [zm2,~]=CreateZoneMat(size(zm2), PF_rad, max_inds2, orig_rates2);
        
        rates_2= nan(1,length(max_inds1));
        for cen= 1:length(max_inds1);
            rates_2(cen)= zm2(max_inds1(cen,1), max_inds1(cen,2));
        end
        
        corr_firing= corrcoef(rates_1,rates_2);
        all_corrs(i)=corr_firing(2);
        
        pr1_orig=rates_1;
        pr2_orig=rates_2;
        
        rates_1(pr1_orig==0 | pr2_orig==0)= [];
        rates_2(pr1_orig==0 | pr2_orig==0)= [];
        
        if length(rates_1) >2
            corr_two= corrcoef(rates_1,rates_2);
            all_corrs2(i)=corr_two(2);
        else
            all_corrs2(i)=nan;
        end
    end
    
    mean_corr(k)=nanmean(all_corrs);
    mean_corr2(k)=nanmean(all_corrs2);
    
    
    save('shuffled same areans firing stability analysis', ...
    'mean_corr','mean_corr2');
    k
    
end


% actual mean corr= 0.5838
% actual mean corr2= 0.5538
