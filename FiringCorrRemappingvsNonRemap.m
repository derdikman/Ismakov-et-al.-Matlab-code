function FiringCorrRemappingvsNonRemap

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis');
load('remapping data info.mat')

% find arenas that are the same context
same_arena_inds=cell(1,length(arena_types_all)); 
for h=1:length(arena_types_all)
    types=arena_types_all{h};
    for a1=1:length(types)-1
        for a2=a1+1:length(types)
            type1=types{a1};
            type2=types{a2};
            if type1(1:2)==type2(1:2)
        same_arena_inds{h}=[a1,a2];
            end
        end
    end
end


PF_thresh=3; % at least 3 fields
corr_thresh= 0.3;
phase_thresh= 0.2;

remap_count=1;
no_remap_count=1;
all_count=1;


for h= 1:length(zone_mats_all)

    % trial= inds(h);
    trial=h;

    % open all arena parameters
    %     hyper_inds_all_arenas= norm_max_index_all{trial};
    max_inds_all_arenas= max_indices_all{trial};
    rates_all_arenas= peak_rates_all{trial};
    zone_mats_all_arenas=zone_mats_all{trial};
    stability_inds= same_arena_inds{trial};
    rate_maps_all_arenas= rate_mats_all{trial};
    MD_scores_all_arenas= MD_scores_all{trial};
    gridness_all_arenas=gridness_all{trial};
    autocorrs_all_arenas= autocorrs_all{trial};

    % for rescaling arenas, remove below:
    % check if same context trials is stable enough
    rate_map_1= rate_maps_all_arenas{stability_inds(1)};
    rate_map_2= rate_maps_all_arenas{stability_inds(2)};
    zone_mat_1= zone_mats_all_arenas{stability_inds(1)};
    zone_mat_2= zone_mats_all_arenas{stability_inds(2)};

    [zone_mat_1, zone_mat_2]= ArenaSameSize(zone_mat_1, zone_mat_2);

    % convert zone mat to 0s and 1s
    zone_mat_1(zone_mat_1~=0)=1;
    zone_mat_2(zone_mat_2~=0)=1;

  
    %xcorr= Cross_Correlation(zone_mat_1,zone_mat_2);
    %R_outer=FindROuter(xcorr);
    %gridness3= GridnessRadiusXcorr(xcorr,R_outer);

    corr_2d=corr2(zone_mat_1,zone_mat_2); 
    if corr_2d > corr_thresh % if same context stable enough
        %
        %         acorr1= Cross_Correlation(rate_map_1,rate_map_1);
        %         spacing1= findPlaceFieldRadius(acorr1) * (10/7);
        %         acorr2= Cross_Correlation(rate_map_2,rate_map_2);
        %         spacing2= findPlaceFieldRadius(acorr2) * (10/7);
        %
        %         [phase_shift]= FindPhaseShift(xcorr, [spacing1 spacing2]);
        %
        %         if phase_shift < phase_thresh % if same context trial stable enough

        len=length(rate_maps_all_arenas);

        for a= 1:len-1
            for b=a+1:len

                if ~isequal([a b], stability_inds)

                rate_mat_1= rate_maps_all_arenas{a};
                rate_mat_2= rate_maps_all_arenas{b};

                max_inds_1= max_inds_all_arenas{a};
                max_inds_2= max_inds_all_arenas{b};

                rates_1=rates_all_arenas{a};
                rates_2=rates_all_arenas{b};

                zone_mat_1= zone_mats_all_arenas{a};
                zone_mat_2= zone_mats_all_arenas{b};

                %                 hyper_ind_1= hyper_inds_all_arenas(a,:);
                %                 hyper_ind_2= hyper_inds_all_arenas(b,:);

                [lenn_1,~] =size(max_inds_1);
                [lenn_2,~]= size(max_inds_2);

                %normalize max inds to arena size
                %                 size_2= size(rate_mat_2);
                %                 norm_max_inds_2=nan(lenn_2,2);
                %                 norm_max_inds_2(:,1)= max_inds_2(:,1)/size_2(1);
                %                 norm_max_inds_2(:,2)= max_inds_2(:,2)/size_2(2);



                %for remapping:
                MD_1= MD_scores_all_arenas(a);
                MD_2= MD_scores_all_arenas(b);

                grid1=gridness_all_arenas(a);
                grid2=gridness_all_arenas(b);

                acorr1= autocorrs_all_arenas{a};
                acorr2= autocorrs_all_arenas{b};

                % if enough number of fields in both trials and all other
                % parameters

                if grid1>=0.3 || grid2>=0.3

                    if  MD_1 <= 0.25 || MD_2 <= 0.25

                        if lenn_1 >= PF_thresh && lenn_2 >= PF_thresh

                            [zone_mat_1, zone_mat_2]= ...
                                ArenaSameSize(zone_mat_1, zone_mat_2);

                            % need only for images:
                            %zone_mat_1_orig= zone_mat_1;
                            %zone_mat_2_orig= zone_mat_2;

                            % convert zone mat to matrix of 0s and 1s
                            zone_mat_1(zone_mat_1~=0)=1;
                            zone_mat_2(zone_mat_2~=0)=1;

                            % find phase shift
                            xcorr= Cross_Correlation(rate_mat_1,rate_mat_2);

                            spacing1= findPlaceFieldRadius(acorr1) * (10/7);
                            spacing2= findPlaceFieldRadius(acorr2) * (10/7);

                            [phase_shift]= FindPhaseShift(xcorr, [spacing1 spacing2]);

                            %xcorr= Cross_Correlation(zone_mat_1,zone_mat_2);
                            %R_outer=FindROuter(xcorr);
                            %[gridness3] = GridnessRadiusXcorr(xcorr,R_outer);

                             corr_2d=corr(zone_mat_1,zone_mat_2);
                             corr_2d=corr_2d(2);
                            [gc_corr, orig, next]= FindFiringProfileGC(max_inds_1, max_inds_2, rates_1, rates_2);

                            all_corrcoefs(all_count)=gc_corr;
                            phase_shifts(all_count)= phase_shift;
                            %xcorr_scores(all_count)=gridness3;
                            corr_scores(all_count)=corr_2d;
                            
                            rms{1}=rate_mat_1;
                            rms{2}=rate_mat_2;
                            all_rms{all_count}=rms;
                                zms{1}=zone_mat_1;
                                zms{2}=zone_mat_2;
                                all_zms{all_count}=zms;
                            all_count=all_count+1;

                            %divide by category

                            if phase_shift < phase_thresh & corr_2d > corr_thresh %didnt remap

                                no_remap_corrs(no_remap_count)=gc_corr;
                                no_remap_count=no_remap_count+1;
                                
                                
                            elseif phase_shift >= phase_thresh | corr_2d <= corr_thresh

                                remap_corrs(remap_count)=gc_corr;
                                remap_count=remap_count+1;

                            end


                            disp('');
                        end
                        end %of if same context
                    end %of if both arenas pass all criteria
                end % of arena 2
            end %of arena 1
            disp('');
        end %end of if same context phase stable
        %  end % end of if same context xcorr gridness stable
    end % of all cells
 end

    save('remapping vs noremapping corrs corr threshold 3', 'remap_corrs', 'no_remap_corrs',...
        'corr_scores','phase_shifts','all_corrcoefs','same_arena_inds',...
        'all_rms','all_zms')
