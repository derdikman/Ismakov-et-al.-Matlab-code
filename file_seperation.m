function file_seperation

parms.dir_load_data = '\\192.114.21.198\Dori_Data\data\rebekkah\data sets\original data sets\Bonnevie 2013 original data';
parms.dir_save_data='\\192.114.21.198\Dori_Data\data\rebekkah\data sets\3 datasets G3MD15PF3';

SeperateFiles(parms,1)

parms.dir_load_data = '\\192.114.21.198\Dori_Data\data\rebekkah\data sets\original data sets\Derdikman Data';

SeperateFiles(parms,2)

parms.dir_load_data = '\\192.114.21.198\Dori_Data\data\rebekkah\data sets\original data sets\sargolini with histology';

SeperateFiles(parms,3)

function SeperateFiles(parms,num)

cd(parms.dir_load_data)
dir_name= parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

grid_thresh=0.3; %0.3 --- 0.2
MD_thresh=0.15;  %0.25 --- 0.3
fields_thresh=3;

% divide into three folders depending on gridness score

%parms.beg_cycle=pi/2;% max point of theta(+0),min ponit (+pi),
%parms.num_of_direction_bins=120;
parms.bin_size=3;
parms.sigma = 1.5;

not_circular_sum=0;
gridness_sum=0;
MD_sum=0;
sum_fields=0;
led_sum=0;

count=0;
count2=0;

for i=1:length(file_names)
    file_name = file_names{i};
    dat=load(file_name);
 
    
    i
    
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
    
  %  if ~isempty(Cell.pos.x2)
        
        led_sum=led_sum+1;
        
%        % find the mean position of the two LEDs
    if isfield(Cell.pos,'x2')&& ~isempty(Cell.pos.x2) 
     pos_mean_x=(Cell.pos.x+ Cell.pos.x2)/2;
     pos_mean_y=(Cell.pos.y+ Cell.pos.y2)/2;
     
     count2=count2+1;
    else
       pos_mean_x=(Cell.pos.x); % + Cell.pos.x2)/2;
       pos_mean_y=(Cell.pos.y); % + Cell.pos.y2)/2;
       
       count=count+1;
    end

        % build the axis of location when spikes are made
        spk_x=interp1(Cell.pos.t,pos_mean_x,Cell.spk.t);
        spk_y=interp1(Cell.pos.t,pos_mean_y,Cell.spk.t);
        
        % Create Rate Maps (smoothed and unsmoothed)
        rate_mat= CreateRateMap(pos_mean_x,pos_mean_y,Cell.pos.t,spk_x,spk_y,Cell.spk.t,parms);
        
        autocorr = Cross_Correlation(rate_mat,rate_mat);
        R_outer = FindROuter(autocorr);
        gridness = GridnessRadius(autocorr,R_outer,i);
        
        if gridness >= grid_thresh
            
            gridness_sum= gridness_sum+1;
            
            % calculate HD
            %          [~,HD_score,~]=ComputeHeadDirectionality...
            %              (Cell.pos.t,Cell.pos.x, Cell.pos.y,Cell.pos.x2,Cell.pos.y2,Cell.spk.t,parms);
            
            [~,MD_score,~]=ComputeMovingDirectionalityWithSpeedThreshHold...
                (Cell.pos.t,Cell.pos.x,Cell.pos.y,Cell.spk.t,2);
            
            if MD_score <= MD_thresh
                
                MD_sum=MD_sum+1;
                
                max_inds=FindMaxIndsRateMap(rate_mat);
                
                auto_max_inds=FindAutoMaxInds(autocorr);
                PF_radius=findPlaceFieldRadius(autocorr, auto_max_inds);
                max_inds=RemoveTooCloseMaxInds(max_inds, PF_radius, rate_mat, 1.9);
                
                [number_of_fields,~]= size(max_inds);
                
                if number_of_fields >= fields_thresh
                    
                    sum_fields=sum_fields+1;
                    
                    rm_size=size(rate_mat);
                    percent_nan= sum(sum(isnan(rate_mat)))/(rm_size(1)*rm_size(2));
                    if percent_nan < 0.13 % if not a circular arena
                        
                        not_circular_sum=not_circular_sum+1;
                        
                        rm_sizes(not_circular_sum,:)=rm_size;
                        
                        copyfile(file_name, parms.dir_save_data);
                        cd(parms.dir_load_data);
                        
                    end
                end % if enough fields
            end % if non-HD
        end % if grid-y enough
   % end % if non-circular
end

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info')
save(sprintf('files that passed criteria 3 sets %d', num), 'sum_fields', 'MD_sum', ...
    'not_circular_sum', 'gridness_sum', 'rm_sizes','led_sum')
