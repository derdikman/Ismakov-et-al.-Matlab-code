function PF_radius= findPlaceFieldRadius(autocorr, auto_max_inds)

if ~exist('auto_max_inds')
auto_max_inds= FindAutoMaxInds(autocorr);
end

[size_x, size_y]= size(autocorr);

 [len, ~] = size(auto_max_inds);
    
        cen = 1:len;
         distances = Distance(auto_max_inds(cen, 1),auto_max_inds(cen, 2),(size_x/2)+0.5,(size_y/2)+0.5);
        
        new_distances = sort(distances (:));

        PF_radius = new_distances(2)/2; %% finds 2nd min since min is 0

        PF_radius= 0.7*PF_radius;