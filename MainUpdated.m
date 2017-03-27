function MainUpdated
%Date: 19 of Jan, 2016. Updated by Rebekkah.

dbstop if error

parms.dir_load_data = '\\192.114.21.198\Dori_Data\data\rebekkah\data sets\3 datasets G3MD15PF3';
%parms.dir_save_pictures= 'N:\users\rebekkah\results and info of analysis\final images';

cd(parms.dir_load_data);

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

len= length(file_names); 

fanos=nan(1,len);
cvs=nan(1,len);
norm_dist=nan(1,len);
norm_second_dist=nan(1,len);
rm_size=nan(len,2);
number_zone_mat= cell(1,len);
max_over_means=nan(1,len);
rate_mats_all= cell(1,len);
zone_mats_all= cell(1,len);
autocorrs_all=cell(1,len);
histology=cell(1,len);
PF_radii=nan(1,len);
max_indices=cell(1,len);
norm_max_index=nan(len,2);
peak_rates_all=cell(1,len);

pos_x_all=cell(1,len);
pos_y_all=cell(1,len);
pos_t_all=cell(1,len);
spk_x_all=cell(1,len);
spk_y_all=cell(1,len);
spk_t_all=cell(1,len);

% enumerate on cells
for i =1:len
    file_name = file_names{i};
    dat = load(file_name);
    
    if ~isfield(dat,'S')
        %for Bonnevie data:
        Cell=dat.db.B(1);
        Cell.pos = Cell.pos_data;
        Cell.pos.x=Cell.pos.x1;
        Cell.pos.y=Cell.pos.y1;
        Cell.spk= Cell.spike_data;
        Cell.spk.t= Cell.spk.ts;
    else
        Cell=dat.S;
    end
    
    % find the mean position of the two LEDs
    if isfield(Cell.pos,'x2')&& ~isempty(Cell.pos.x2)
        pos_mean_x=(Cell.pos.x+ Cell.pos.x2)/2;
        pos_mean_y=(Cell.pos.y+ Cell.pos.y2)/2;
    else
        pos_mean_x=(Cell.pos.x); % + Cell.pos.x2)/2;
        pos_mean_y=(Cell.pos.y); % + Cell.pos.y2)/2;
    end
    
    
    % build the axis of location when spikes are made
    spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
    spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
    
    % Create Rate Maps (smoothed and unsmoothed)
    parms.sigma=1.5;
    parms.bin_size=3;
    rate_mat= CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
    
    %parms.bin_size=6;
    %rate_mat_ns= CreateRateMapNoSmooth(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
    
    % Create AutoCorr
    autocorr=Cross_Correlation(rate_mat, rate_mat);
    
    %Find AutoMaxInds
    auto_max_inds= FindAutoMaxInds(autocorr);
    
    % Find PF_radius
    PF_radius=findPlaceFieldRadius(autocorr, auto_max_inds);
    
    % Find Max_Inds of smoothed & nonsmoothed rate mat
    max_inds= FindMaxIndsRateMap(rate_mat);
    
    % plot max inds before filtering:
    if i==3 || i==4 || i==5
        im=[1,3,5];
    subplot(3,2,im(i-2));
    imagesc(rate_mat);axis xy; axis off; hold on;
    plot(max_inds(:,2),max_inds(:,1), 'kx','Markersize',5, 'linewidth',5)
    end
    
    strength= 1.9;
    max_inds= RemoveTooCloseMaxInds(max_inds, PF_radius, rate_mat, strength);
    
     if i==3 || i==4 || i==5
     subplot(3,2,im(i-2)+1);
    imagesc(rate_mat);axis xy;axis off; hold on;
    plot(max_inds(:,2),max_inds(:,1), 'kx','Markersize',5,'linewidth',5)
     end
    
%     PF_radius_ns= PF_radius/2; %half since bin_size is 6 instead of 3
%     max_inds_ns= FindMaxIndsRateMap(rate_mat_ns);
%     max_inds_ns= RemoveTooCloseMaxInds(max_inds_ns, PF_radius_ns, rate_mat_ns, strength);
    
    % Find peak firing rate at max_inds
    peak_rates=findPeakRates(max_inds,rate_mat);
    
    % Find Fano factor
    fanos(i)= var(peak_rates)/mean(peak_rates);
    cvs(i)=std(peak_rates)/mean(peak_rates);
    
    % Find location of max-firing field (hyperfield)
    ind= find(peak_rates==max(peak_rates));
    ind=ind(1); %in cases where same rate, take first
    max_index= max_inds(ind,:);
    size_rate_mat= size(rate_mat);
    norm_max_index(i,:)= max_index ./ size_rate_mat;
    
    %max_indices_ns{i}= max_inds_ns; %save for future max peak shuffling
    max_indices{i}= max_inds; %save for future max peak shuffling
    
    %find distance of max field location to border
    [~,norm_dist(i)]= findDistPtToBorder([1 size_rate_mat(1)], [1 size_rate_mat(2)], max_index);
    
    %save rate mat sizes and nonsmooth peak rates for future max peak shuffling
    rm_size(i,:)= size_rate_mat;
   % peak_rates_ns_all{i}=peak_rates_ns;
    
    peak_rates_all{i}=peak_rates; %save for border vs nonborder mean rates analysis
    
    size_rate_mat= size(rate_mat);
    %     %find distance of second max field location to border
    peak_rates_wo_max= peak_rates;
    peak_rates_wo_max(peak_rates_wo_max==max(peak_rates_wo_max))= nan;
    second_ind= find(peak_rates_wo_max==max(peak_rates_wo_max));
    second_ind=second_ind(1); %in case same rate take first
    second_max_index= max_inds(second_ind,:);
    [~,norm_second_dist(i)]= findDistPtToBorder([ 1 size_rate_mat(1)], [1 size_rate_mat(2)], second_max_index);
    
    % create zone mats for images and number_zone_mat for future shuffling
    [zone_mat, number_zone_mat{i}]= CreateZoneMat(size(rate_mat), PF_radius, max_inds, peak_rates);
    %[zone_mat_ns, number_zone_mat{i}]= CreateZoneMat(rate_mat_ns, PF_radius_ns, max_inds_ns, peak_rates_ns);
    
    sr=sort(peak_rates);
    max_over_means(i)= sr(end)/ mean(sr(1:end-1));
    
    %image
    %     fig=figure;
    %     n=2;m=4;
    %     subplot(n,m,1)
    %     plot(pos_mean_x,pos_mean_y,'k');hold on;
    %     plot(spk_x,spk_y,'.r');
    %     axis equal;axis off;
    %     axis ij;
    %     title_name= file_name(14:end-4);
    %     title(sprintf('%s', title_name), 'Interpreter', 'none');
    %
    %     subplot(n,m,2)
    %     imagesc(rate_mat);
    %     axis equal; axis off;
    %     title(sprintf('%0.1f Hz', max(peak_rates)), 'HorizontalAlignment', 'left');
    %
    %     subplot(n,m,3)
    %     imagesc(zone_mat);
    %     axis equal; axis off;
    %
    %     subplot(n,m,4)
    %     plot(1:length(peak_rates),sort(peak_rates), 'ko-');
    %     title(sprintf('Fano factor=%0.1f', fano_factor(i)));
    %
    %     subplot(n,m,5)
    %     imagesc(autocorr);
    %     axis equal; axis off;
    %
    %     subplot(n,m,6)
    %     imagesc(rate_mat_ns);
    %     axis equal; axis off;
    %     title(sprintf('%0.1f Hz', max(peak_rates_ns)),'HorizontalAlignment', 'left');
    %
    %     subplot(n,m,7)
    %     imagesc(zone_mat_ns);
    %     axis equal; axis off;
    %     title(sprintf('dist=%0.2f', norm_dist(i)));
    %
    %     subplot(n,m,8)
    %     plot(1:length(peak_rates_ns),sort(peak_rates_ns), 'ko-');
    %     title(sprintf('Fano factor=%0.1f', fano_factor_ns));
    %
    %     %save images
    %     cd(parms.dir_save_pictures);
    %     saveas(fig,sprintf('%s.jpg',file_name)); %
    %     cd(parms.dir_load_data);
    
    rate_mats_all{i}= rate_mat;
    %rate_mats_ns_all{i}= rate_mat_ns;
    zone_mats_all{i}= zone_mat;
    %zone_mats_ns_all{i}= zone_mat_ns;
    autocorrs_all{i}=autocorr;
    PF_radii(i)=PF_radius;
    
    pos_x_all{i}=pos_mean_x;
    pos_y_all{i}=pos_mean_y;
    pos_t_all{i}=Cell.pos.t;
    spk_x_all{i}=spk_x;
    spk_y_all{i}=spk_y;
    spk_t_all{i}=Cell.spk.t;
    
    if isfield(dat,'db') %Bonnevie
        histology{i}=dat.db.area;
    elseif isfield('histology_layer',Cell)
        histology{i}=Cell.histology_layer;
    else
        histology{i}=nan;
    end
    
    R_outer = FindROuter(autocorr);
    gridness = GridnessRadius(autocorr,R_outer,i);
    gridness_all(i)=gridness;
    [~,MD_score,~]=ComputeMovingDirectionalityWithSpeedThreshHold...
                (Cell.pos.t,Cell.pos.x,Cell.pos.y,Cell.spk.t,2);
     MD_all(i)=MD_score;
    
    
   % close all
end

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis');
% save('3sets G3MD25 data and results', 'fanos', 'norm_dist', ...
%     'norm_max_index', 'max_indices_ns', 'norm_second_dist', ...
%     'rm_size', 'peak_rates_ns_all', 'number_zone_mat', ...
%     'peak_rates_sm_all', 'number_zone_mat_sm', 'max_indices_sm', ...
%     'max_over_means', 'rate_mats_all', 'rate_mats_ns_all', ...
%     'zone_mats_all', 'zone_mats_ns_all', 'autocorrs_all',...
%     'histology', 'PF_radii'...
%     'pos_x_all','pos_y_all','pos_t_all',...
%     'spk_x_all','spk_y_all', 'spk_t_all');

save('3sets G3MD15PF3 data and results again', 'fanos', 'cvs','norm_dist', ...
    'norm_max_index', 'max_indices', 'norm_second_dist', ...
    'rm_size', 'number_zone_mat', ...
    'peak_rates_all', 'file_names',  ...
    'max_over_means', 'rate_mats_all', ...
    'zone_mats_all', 'autocorrs_all',...
    'histology', 'PF_radii',...
    'pos_x_all','pos_y_all','pos_t_all',...
    'spk_x_all','spk_y_all', 'spk_t_all', 'gridness_all','MD_all');

disp('');
