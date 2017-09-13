function [fano_factor, peak_rates,cv,max_inds]= findFano(rate_mat,PF_radius)

% Find Max_Inds of smoothed & nonsmoothed rate mat
max_inds= FindMaxIndsRateMap(rate_mat);
strength= 1.9;
max_inds= RemoveTooCloseMaxInds(max_inds, PF_radius, rate_mat, strength);

% Find peak firing rate at max_inds
peak_rates = findPeakRates(max_inds, rate_mat);

% Find Fano factor
if length(peak_rates) >1
fano_factor= var(peak_rates)/mean(peak_rates);
cv= std(peak_rates)/mean(peak_rates);
else
fano_factor= nan; 
cv=nan;
end