% Date: 20 of Jan 2016
dbstop if error

parms.dir_load_data = '\\192.114.21.198\Dori_Data\data\rebekkah\data sets\rescaling arena data set kate jefferey';
cd(parms.dir_load_data);

parms.bin_size=5;
parms.sigma=3;

dir_name=parms.dir_load_data;
dir_list=dir(strcat(dir_name,'\*.mat'));
file_names={dir_list.name};

count=1;
norm_hyperinds=[];

number_of_cells=0;
grid_num=0;

for i =1:length(file_names) %open all files
    file_name = file_names{i};
    load(file_name);
    
    gridness_scores= nan(1,length(tint));
    max_inds=cell(1,length(tint));
    rm_size= cell(1,length(tint));
    rate_mats= cell(1,length(tint));
    zone_mats= cell(1,length(tint));
    max_indices= cell(1,length(tint));
    peak_rates_arena= cell(1,length(tint));
    autocorrs= cell(1,length(tint));
    MD_scores_all_arenas=nan(1,length(tint));
    
    pos_x_aa=cell(1,length(tint));
    pos_y_aa=cell(1,length(tint));
    pos_t_aa=cell(1,length(tint));
    spk_x_aa=cell(1,length(tint));
    spk_y_aa=cell(1,length(tint));
    spk_t_aa=cell(1,length(tint));
    
    % find posx and posy
    
    for cell_num=1:length(cells) %open every cell per recording
        
        number_of_cells=number_of_cells+1;
        
        mySpikes=[];
        
        for arena_count= 1:length(tint); %open every arena
            
            if tets(cell_num)==1
                mySpikes=find(cutTet1{arena_count}==cells(cell_num));
            elseif tets(cell_num)==2
                mySpikes=find(cutTet2{arena_count}==cells(cell_num));
            elseif tets(cell_num)==3
                mySpikes=find(cutTet3{arena_count}==cells(cell_num));
            elseif tets(cell_num)==4
                mySpikes=find(cutTet4{arena_count}==cells(cell_num));
            end
            
            if length(tint(arena_count).tetrode)>= tets(cell_num) %why is tetrode 4 missing when spikes are supposedly found at tetrode 4? file_name=216_6_6_12
                
                len= length(tint(arena_count).tetrode(tets(cell_num)).pos_sample); % why is mySpikes longer than the length of pos_sample???
                mySpikes(mySpikes>len)=[]; % why is mySpikes longer than the length of pos_sample???
                mySpikes=mySpikes';
                
                myPosSample= tint(arena_count).tetrode(tets(cell_num)).pos_sample(mySpikes);
                
                spk_x=tint(arena_count).pos.xy(myPosSample,1);
                spk_y=tint(arena_count).pos.xy(myPosSample,2);
 
                pos_t=0.02:0.02:length(tint(arena_count).pos.xy)*0.02;
 
                pos_x= tint(arena_count).pos.xy(:,1);
                pos_y= tint(arena_count).pos.xy(:,2);
                
                spk_t=nan(1,length(spk_x));
                                
                for h=1:length(spk_x)
                   if ~isnan(spk_x(h))
                        spk_xy_ind= find(pos_x==spk_x(h) & pos_y==spk_y(h));
                        spk_t(h)= pos_t(spk_xy_ind(1));
                   end
                end
                
                pos_x_aa{arena_count}=pos_x;
                pos_y_aa{arena_count}=pos_y;
                pos_t_aa{arena_count}=pos_t;
                spk_x_aa{arena_count}=spk_x;
                spk_y_aa{arena_count}=spk_y;
                spk_t_aa{arena_count}=spk_t;
                
                % create rate mat 
                [rate_mat]=CreateRateMap(pos_x,pos_y,pos_t,spk_x,spk_y,spk_t,parms);
                
                % create autocorr and find gridness
                autocorr = Cross_Correlation(rate_mat,rate_mat);
                R_outer = FindROuter(autocorr);
                [gridness] = GridnessRadius(autocorr,R_outer,i);
                
                gridness_scores(arena_count) = gridness;
                
                %save all info for future shuffling
                max_inds= FindMaxIndsRateMap(rate_mat);
                auto_max_inds= FindAutoMaxInds(autocorr);
                PF_radius=findPlaceFieldRadius(autocorr, auto_max_inds);
                max_inds= RemoveTooCloseMaxInds(max_inds, PF_radius, rate_mat, 1.9);
                
                rm_size{arena_count}= size(rate_mat);
                
                rate_mats{arena_count}= rate_mat;
                max_indices{arena_count}= max_inds;
                autocorrs{arena_count}=autocorr;
                
                peak_rates= findPeakRates(max_inds,rate_mat);
                
                [zone_mats{arena_count},~]= CreateZoneMat(size(rate_mat), PF_radius, max_inds, peak_rates);
                
                peak_rates_arena{arena_count}= peak_rates;
                
         
               [~,MD_scores_all_arenas(arena_count),~]=ComputeMovingDirectionalityWithSpeedThreshHold...
                (pos_t,pos_x,pos_y,spk_t,0);
                
            end
        end %end of arena_count
               
        max_inds_use= [];
        norm_inds= cell(1,length(gridness_scores));
        if any(gridness_scores > 0.3) %if gridness score >0.3 in at least 1 arena
            
            grid_num=grid_num+1;
            
            if any(MD_scores_all_arenas<0.15) 
            
            % save info if at least one arena of gridness >= 0.3
            rate_mats_all{count}= rate_mats;
            zone_mats_all{count}= zone_mats;
            max_indices_all{count}= max_indices;
            peak_rates_all{count}= peak_rates_arena;
            autocorrs_all{count}= autocorrs;
            gridness_all{count}= gridness_scores;
            rm_sizes_all{count}=rm_size;
            MD_scores_all{count}=MD_scores_all_arenas;
            
            pos_x_all{count}=pos_x_aa;
            pos_y_all{count}=pos_y_aa;
            pos_t_all{count}=pos_t_aa;
            spk_x_all{count}=spk_x_aa;
            spk_y_all{count}=spk_y_aa;
            spk_t_all{count}=spk_t_aa;
            
            shape_seqs{count}=shapeSeq;
            
            for h= 1:length(max_indices)
                max_inds_use= (max_indices{h});
                size_use= rm_size{h};
                max_inds_use(:,1)=max_inds_use(:,1)/size_use(1);
                max_inds_use(:,2)=max_inds_use(:,2)/size_use(2);
                norm_inds{h}= max_inds_use;
            end
            
            all_norm_inds{count}= norm_inds;
            
            title=sprintf('%s cell %d', file_name, cell_num);
            filenames{count}= title;
            
            count=count+1;
            end  % of if at least one MD score <= 0.25
        end % end of if at least one 0.3 gridness score  
    end %end of cell_num
end %end of filenames

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis');
save('rescaling arenas data info', 'rate_mats_all', 'zone_mats_all', ...
    'max_indices_all', 'peak_rates_all', 'autocorrs_all', 'filenames', ...
    'norm_hyperinds', 'gridness_all', 'rm_sizes_all',...
    'pos_x_all', 'pos_y_all', 'pos_t_all', ...
    'spk_x_all', 'spk_y_all', 'spk_t_all', 'MD_scores_all', 'grid_num',...
    'number_of_cells', 'shape_seqs',...
    'all_norm_inds')