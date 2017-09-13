function [auto_max_inds] = FindAutoMaxInds(autocorr)

[size_x, size_y] = size(autocorr);

auto_max_inds_len = 0;

for fig_i = 2:size_x-1
    for j = 2:size_y-1
        if autocorr(fig_i,j) > autocorr(fig_i+1,j) && ...
                autocorr(fig_i,j) > autocorr(fig_i-1,j) && ...
                autocorr(fig_i,j) > autocorr(fig_i,j+1) && ...
                autocorr(fig_i,j) > autocorr(fig_i,j-1)
            % plot(j,fig_i,'x');
            
            auto_max_inds_len = auto_max_inds_len+1;
            auto_max_inds(auto_max_inds_len,1) = fig_i;   %indices of maximum pts
            auto_max_inds(auto_max_inds_len,2) = j;
        end
    end
end

