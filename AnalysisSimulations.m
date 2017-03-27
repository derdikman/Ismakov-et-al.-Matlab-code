function AnalysisSimulations

% simulates spike trains using rate map and trajectory
% finds all information using simulated rate maps to compare with original

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis')
load('3sets G3MD15PF3 data and results.mat', ...
    'pos_x_all','pos_y_all','pos_t_all', 'zone_mats_all','max_indices')
load('G3MD15PF3 gaussian mats and PF rads larger size.mat')

% parameters for rate map smoothing
parms.bin_size= 3;
parms.sigma= 1.5;

num= 100; % number of simulations per cell

len= length(pos_x_all); % number of cells in set

% initialize all variables
all_fano_factors= nan(len,num);
all_max_over_means= nan(len,num);
all_stability_corrs= nan(len,num);
all_cvs= nan(len,num);

for i= 1:len %go through all cells    
    pos_x= pos_x_all{i};
    pos_y= pos_y_all{i};
    pos_t= pos_t_all{i};
    
    gaussian_mat= gaussian_mats{i};
    PF_rad= PF_radii(i);
    zone_mat= zone_mats_all{i};
    max_inds= max_indices{i};
    
    % simulates spike trains using Gaussian map and original trajectory
    [spk_t]= Simulate_Spike_Train(pos_x,pos_y,pos_t,gaussian_mat,parms,num);
    
    % initialize all variables for cell
    simulated_fano= nan(1,num);
    max_over_mean= nan(1,num);
    all_corrs= nan(1,num);
    simulated_cvs= nan(1,num);
    
    for h= 1:num % goes through every simulation run    
        spk_t_use= spk_t(h,:); 
        spk_t_use(isnan(spk_t_use))= [];
        
        spk_x= interp1(pos_t,pos_x,spk_t_use);
        spk_y= interp1(pos_t,pos_y,spk_t_use);
        
        % create simulated rate map        
        simulated_rate_map= CreateRateMap(pos_x, pos_y, pos_t, ...
            spk_x, spk_y, spk_t_use, parms);
        
        % firing stability of split session
        [~,all_corrs(h), ~,~,~,~,~,~,~,~,~,~]= ...
        FiringStabilitySplitSessions(pos_x,pos_y,pos_t,spk_t_use,...
            zone_mat, max_inds,parms);
                
        [simulated_fano(h), peak_rates,simulated_cvs(h),~]= ...
            findFano(simulated_rate_map,PF_rad);
              
        sr= sort(peak_rates);
        max_over_mean(h)= sr(end)/ mean(sr(1:end-1));
      
    end % end of simulation run
    
    all_fano_factors(i,:)= simulated_fano;
    all_max_over_means(i,:)= max_over_mean;
    all_stability_corrs(i,:)= all_corrs;
    all_cvs(i,:)= simulated_cvs;
    
    save('simulated results G3MD15PF3 100x larger field sizes.mat', ...
       'all_cvs','all_fano_factors','all_max_over_means','all_stability_corrs')
        
end % end of all cell simulations
