function [zone_mat, number_zone_mat]= CreateZoneMat(size_rm, PF_radius, max_inds, peak_rates)

% create zone mat with peak rates at field centers
zone_mat= zeros(size_rm);
number_zone_mat= zeros(size_rm);

[max_inds_len,~]= size(max_inds);

for cen=1:max_inds_len
    zone_mat(max_inds(cen,1), max_inds(cen,2))= peak_rates(cen);
end

% if distance to max pt is less than PF radius, assign it value at max pt

for cen=1:max_inds_len;
    for fig_i =1:size_rm(1)
        for j =1:size_rm(2)
            if Distance(fig_i, j, max_inds(cen,1), max_inds(cen,2)) < PF_radius  %change this depending on how large you want fields to be
                zone_mat(fig_i,j)= peak_rates(cen);
                number_zone_mat(fig_i,j)= cen;
            end
        end
    end
end
