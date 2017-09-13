function [max_inds] = RemoveTooCloseMaxInds(max_inds, PF_radius, rate_mat, strength)

h=1;
too_close=[];
remove=[];

[max_inds_len,~]= size(max_inds);

for cen = 1:max_inds_len-1
    for cen2= (cen+1):max_inds_len
        peak_distance = Distance(max_inds(cen,1), max_inds(cen,2),max_inds(cen2,1), max_inds(cen2,2));
        if peak_distance < strength *(PF_radius)      %determine what is too close, 1.25 too small
            too_close (h,1) = cen;
            too_close(h,2) = cen2;
            
            h=h+1;
        end
    end
end

if ~isempty (too_close)
    
    [too_close_len, ~] = size(too_close);
    
    %remove one of the peaks that are too close with lower firing rate
    
    h=1;
    
    for cen = 1: too_close_len
        if rate_mat((max_inds((too_close(cen,1)), 1)), (max_inds((too_close(cen,1)), 2))) >=...
                rate_mat((max_inds((too_close(cen,2)), 1)), (max_inds((too_close(cen,2)), 2)))  ;
            
            remove(h) = too_close(cen,2);
            
            h= h+1;
            
        elseif rate_mat((max_inds((too_close(cen,1)), 1)), (max_inds((too_close(cen,1)), 2))) <...
                rate_mat((max_inds((too_close(cen,2)), 1)), (max_inds((too_close(cen,2)), 2)))  ;
            
            remove(h)= too_close(cen,1);
            
            h=h+1;
            
        else
            disp ('wtf')
        end
    end
    
    max_inds(remove,:) = [];

end
