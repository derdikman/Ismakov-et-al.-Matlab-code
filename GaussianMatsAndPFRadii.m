function GaussianMatsAndPFRadii

load('3sets G3MD15PF3 data and results.mat')

gaussian_mats=cell(1,length(fanos));
for i =1:length(fanos)
    
    pos_x=pos_x_all{i};
    pos_y=pos_y_all{i};
    pos_t=pos_t_all{i};
    spk_x=spk_x_all{i};
    spk_y=spk_y_all{i};
    spk_t=spk_t_all{i};
    
    max_inds=max_indices{i};
    rate_mat=rate_mats_all{i}; 
    peak_rates=peak_rates_all{i};
    
   
    gaussian_mats{i}= createGaussianMat(rate_mat, max_inds, mean(peak_rates));

 disp('')
end

save('G3MD15PF3 gaussian mats and PF rads larger size', 'gaussian_mats','PF_radii')