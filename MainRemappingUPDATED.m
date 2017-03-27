function MainRemappingUPDATED
%Date: 05 of Apr, 2016. Updated by Rebekkah.

dbstop if error

parms.dir_load_data = '\\192.114.21.198\Dori_Data\data\rebekkah\All rats\all cells with no pos data fixing';
parms.dir_save_pictures= '\\192.114.21.198\Dori_Data\data\rebekkah\All rats\all cells images';
parms.dir_save_images= '\\192.114.21.198\Dori_Data\data\rebekkah\All rats\images original analysis';

cd(parms.dir_load_data);

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

add=0; % for figure subplot order

count=1;

% enumerate on cells
for i =1:375 %length(file_names)
    file_name = file_names{i};
    dat = load(file_name);
    Cell= dat.db;
    
    len=5;
    
    pos_x= (Cell.pos_data.x1); %+Cell.pos_data.x2)/2;
    pos_y= (Cell.pos_data.y1); %+Cell.pos_data.y2)/2;  
    pos_t= Cell.pos_data.t;
    spk_t= Cell.spike_data.ts;
    
    spk_x=interp1(pos_t,pos_x,spk_t);
    spk_y=interp1(pos_t,pos_y,spk_t);
    
    [~,MD_score,~]=ComputeMovingDirectionalityWithSpeedThreshHold...
        (pos_t,pos_x,pos_y,spk_t,0);
    
    % Create Rate Maps (smoothed and unsmoothed)
    parms.sigma=5;
    parms.bin_size=4;
    rate_mat= CreateRateMap(pos_x,pos_y,pos_t,spk_x,spk_y,spk_t,parms);
    
    %     parms.bin_size=6;
    %     rate_mat_ns= CreateRateMapNoSmooth(pos_mean_x,pos_mean_y,pos_t,spk_x,spk_y,spk_t,parms);
    
    % Create AutoCorr
    autocorr=Cross_Correlation(rate_mat, rate_mat);
    
    %Find AutoMaxInds
    auto_max_inds= FindAutoMaxInds(autocorr);
    
    % Find PF_radius
    PF_radius=findPlaceFieldRadius(autocorr, auto_max_inds);
    
    % Find Max_Inds of smoothed & nonsmoothed rate mat
    max_inds= FindMaxIndsRateMap(rate_mat);
    strength= 1.9;
    max_inds= RemoveTooCloseMaxInds(max_inds, PF_radius, rate_mat, strength);
    
    % Find peak firing rate at max_inds
    peak_rates= findPeakRates(max_inds,rate_mat);
    
    % Find location of max-firing field (hyperfield)- nonsmoothed data
    ind= find(peak_rates==max(peak_rates));
    ind=ind(1); %in cases where same rate, take first
    max_index= max_inds(ind,:);
    size_rate_mat= size(rate_mat);
    norm_max_index= max_index ./ size_rate_mat;
    
    % max_indices{i}= max_inds_ns; %save for future max peak shuffling
    max_indices= max_inds; %save for future max peak shuffling
    
    %find distance of max field location to border
    % [~,norm_dist(i)]= findDistPtToBorder([1 size_rate_mat(1)], [1 size_rate_mat(2)], max_index);
    
    %save rate mat sizes and nonsmooth peak rates for future max peak shuffling
    % rm_size(i,:)= size_rate_mat;
    % peak_rates_all{i}=peak_rates_ns;
    
    %peak_rates_all{i}=peak_rates; %save for border vs nonborder mean rates analysis
    
    % size_rate_mat= size(rate_mat);
    %     %find distance of second max field location to border
    %     peak_rates_wo_max= peak_rates;
    %     peak_rates_wo_max(peak_rates_wo_max==max(peak_rates_wo_max))= nan;
    %     second_ind= find(peak_rates_wo_max==max(peak_rates_wo_max));
    %     second_ind=second_ind(1); %in case same rate take first
    %     second_max_index= max_inds(second_ind,:);
    %     [~,norm_second_dist(i)]= findDistPtToBorder([1 size_rate_mat(1)], [1 size_rate_mat(2)], second_max_index);
    
    % create zone mats for images and number_zone_mat for future shuffling
    [zone_mat, ~]= CreateZoneMat(size(rate_mat), PF_radius, max_inds, peak_rates);
    % [zone_mat_ns, number_zone_mat{i}]= CreateZoneMat(rate_mat_ns, PF_radius_ns, max_inds_ns, peak_rates_ns);
    
    %find gridness score
    R_outer = FindROuter(autocorr);
    gridness = GridnessRadius(autocorr,R_outer,i);
    
    add=add+1;
    
    %         image
    %                     if add==1
    %                         fig=figure;
    %                     end
    %
    %                     n=len;m=5;
    %                     subplot(n,m,1+(5*(add-1)))
    %                     plot(pos_mean_x,pos_mean_y,'k');hold on;
    %                     plot(spk_x,spk_y,'.r');
    %                     axis equal;axis off;
    %                     axis ij;
    %                     title_name= Cell.cut_file(1:end-4);
    %                     title(sprintf('%s %s %s %s', title_name, Cell.arena, Cell.tetrode, Cell.cell), 'Interpreter', 'none');
    %
    %                     subplot(n,m,2+(5*(add-1)))
    %                     imagesc(rate_mat);
    %                     axis equal; axis off;
    %                     title(sprintf('%0.1f Hz', max(peak_rates)), 'HorizontalAlignment', 'left');
    %
    %                     subplot(n,m,3+(5*(add-1)))
    %                     imagesc(zone_mat);
    %                     title(sprintf('cell = %d', count));
    %                     axis equal; axis off;
    %
    %                     subplot(n,m,4+(5*(add-1)))
    %                     plot(1:length(peak_rates),sort(peak_rates), 'ko-');
    %                     title(sprintf('Fano factor=%0.1f', fano_factor(i)));
    %                     axis square;
    %
    %                     subplot(n,m,5+(5*(add-1)))
    %                     imagesc(autocorr);
    %                     axis equal; axis off;
    %                     title(sprintf('gridness= %0.2f', gridness));
    %
    %                     subplot(n,m,6)
    %                     imagesc(rate_mat_ns);
    %                     axis equal; axis off;
    %                     title(sprintf('%0.1f Hz', max(peak_rates_ns)),'HorizontalAlignment', 'left');
    %
    %                     subplot(n,m,7)
    %                     imagesc(zone_mat_ns);
    %                     axis equal; axis off;
    %                     title(sprintf('dist=%0.2f', norm_dist(i)));
    %
    %                     subplot(n,m,8)
    %                     plot(1:length(peak_rates_ns),sort(peak_rates_ns), 'ko-');
    %                     title(sprintf('Fano factor=%0.1f', fano_factor_ns));
    
    
    
    %save images
    
    zone_mats{add}= zone_mat;
    arena_types{add}= Cell.arena;
    max_inds_all(add,:)= norm_max_index;
    all_max_indices{add}= max_indices;
    peak_rates_all_arenas{add}=peak_rates;
    rate_mats_all_arenas{add}=rate_mat;
    gridness_scores_all_arenas(add)= gridness;
    MD_scores_all_arenas(add)= MD_score;
    autocorrs_all_arenas{add}= autocorr;
    
    posx_aa{add}= pos_x;
    posy_aa{add}= pos_y;
    post_aa{add}= pos_t;
    spkx_aa{add}= spk_x;
    spky_aa{add}= spk_y;
    spkt_aa{add}= spk_t;
    
    if add==len
        
        %             set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 17.4], 'PaperUnits', 'centimeters', 'PaperSize', [17.4, 17.4])
        %
        %             cd(parms.dir_save_pictures);
        %             saveas(fig,sprintf('%s.jpg', file_name(1:end-4))); %
        %   cd(parms.dir_load_data);
        
        if any(gridness_scores_all_arenas>0.3) && any(MD_scores_all_arenas<0.15)
            
            filenames{count}=file_name;
            zone_mats_all{count}= zone_mats;
            arena_types_all{count}= arena_types;
            norm_max_index_all{count}= max_inds_all;
            max_indices_all{count}= all_max_indices;
            peak_rates_all{count}= peak_rates_all_arenas;
            rate_mats_all{count}= rate_mats_all_arenas;
            gridness_all{count}=gridness_scores_all_arenas;
            MD_scores_all{count}=MD_scores_all_arenas;
            autocorrs_all{count}=autocorrs_all_arenas;
            
            pos_x_all{count}= posx_aa;
            pos_y_all{count}= posy_aa;
            pos_t_all{count}= post_aa;
            spk_x_all{count}= spkx_aa;
            spk_y_all{count}= spky_aa;
            spk_t_all{count}= spkt_aa;
            
            count=count+1;
            
        end
        
        zone_mats= cell(1,len);
        arena_types= cell(1,len);
        max_inds_all= nan(len,2);
        all_max_indices=cell(1,len);
        rate_mats_all_arenas=cell(1,len);
        gridness_scores_all_arenas= nan(1,len);
        MD_scores_all_arenas= nan(1,len);
        autocorrs_all_arenas=cell(1,len);
        
        posx_aa= cell(1,len);
        posy_aa= cell(1,len);
        post_aa= cell(1,len);
        spkx_aa= cell(1,len);
        spky_aa= cell(1,len);
        spkt_aa= cell(1,len);
        
        add=0;
    end
    
    i
    
end

disp('');

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis');
save('remapping data info', 'rate_mats_all', 'zone_mats_all', ...
    'max_indices_all', 'peak_rates_all', 'autocorrs_all', 'filenames', ...
    'gridness_all','MD_scores_all', 'arena_types_all',...
    'pos_x_all', 'pos_y_all', 'pos_t_all', ...
    'spk_x_all', 'spk_y_all', 'spk_t_all')

