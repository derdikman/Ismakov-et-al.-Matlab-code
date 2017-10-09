function [zone_mat, peak_rates, max_inds]= getZones(rate_mat)

autocorr = Cross_Correlation(rate_mat,rate_mat);
max_inds= FindMaxIndsRateMap(rate_mat);
auto_max_inds= FindAutoMaxInds(autocorr);
PF_radius=findPlaceFieldRadius(autocorr, auto_max_inds);
max_inds= RemoveTooCloseMaxInds(max_inds, PF_radius, rate_mat, 1.9);

[len,~]=size(max_inds);
peak_rates= nan(1,len);
for cen= 1:len;
    peak_rates(cen)= rate_mat(max_inds(cen,1), max_inds(cen,2));
end

rm_size=size(rate_mat);
[zone_mat, ~]= CreateZoneMat(rm_size, PF_radius, max_inds, peak_rates);