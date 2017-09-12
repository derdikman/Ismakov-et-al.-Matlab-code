function peak_rates = findPeakRates(max_inds, rate_mat)

[len,~]=size(max_inds);
peak_rates = nan(1,len);
for cen = 1:len
    peak_rates(cen)= rate_mat(max_inds(cen,1), max_inds(cen,2));
end

