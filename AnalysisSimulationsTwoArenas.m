function AnalysisSimulationsTwoArenas

% simulates spike trains using rate maps and trajectory for same arenas
% finds all information using simulated rate maps to compare with original

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis')
load('corr coef results of same arenas COMBINED.mat',...
    'pos_x_all1','pos_y_all1','pos_t_all1', ...
    'pos_x_all2','pos_y_all2','pos_t_all2', ...
    'gaussian_mat1','gaussian_mat2','PF_radii1','PF_radii2')

% parameters for rate map smoothing:
parms.bin_size= 5;
parms.sigma= 3;

num= 200; % number of simulations

len= length(pos_x_all1); % number of cells

% initialize variable
all_stability_corrs= nan(len,num);

for i= 1:len % iterates through every cell    
    pos_x1= pos_x_all1{i};
    pos_y1= pos_y_all1{i};
    pos_t1= pos_t_all1{i};
    pos_x2= pos_x_all2{i};
    pos_y2= pos_y_all2{i};
    pos_t2= pos_t_all2{i};
    
    PF_rad1= PF_radii1(i);
    gm1= gaussian_mat1{i};
    PF_rad2= PF_radii2(i);
    gm2= gaussian_mat2{i};
     
    % simulates spike trains using Gaussian map and original trajectory
    [spk_t1]= Simulate_Spike_Train(pos_x1,pos_y1,pos_t1,gm1,parms,num);
    [spk_t2]= Simulate_Spike_Train(pos_x2,pos_y2,pos_t2,gm2,parms,num);
    
    all_corrs=nan(1,num); % initialize variable
    
    for h= 1:num % iterates through every simulation run       
        spk_t_use1= spk_t1(h,:); 
        spk_t_use1(isnan(spk_t_use1))=[];
        spk_t_use2= spk_t2(h,:);
        spk_t_use2(isnan(spk_t_use2))=[];
        
        spk_x1=interp1(pos_t1,pos_x1,spk_t_use1);
        spk_y1=interp1(pos_t1,pos_y1,spk_t_use1);
        spk_x2=interp1(pos_t2,pos_x2,spk_t_use2);
        spk_y2=interp1(pos_t2,pos_y2,spk_t_use2);
        
        % create simulated rate maps
        simulated_rm1= CreateRateMap(pos_x1, pos_y1, pos_t1, ...
            spk_x1, spk_y1, spk_t_use1, parms);
        simulated_rm2= CreateRateMap(pos_x2, pos_y2, pos_t2, ...
            spk_x2, spk_y2, spk_t_use2, parms);
        
        % firing stability of split session:
        [all_corrs(h)]= TwoArenaFiringStability(...
            simulated_rm1,simulated_rm2,PF_rad1,PF_rad2);
                
    end % end of simulation run
end

% results from all simulations and all cells:
all_stability_corrs(i,:)= all_corrs; 

save('simulated results same arenas 200x per cell dori.mat', ...
    'all_stability_corrs')
