function [dist, norm_dist]= findDistPtToBorder(size_x, size_y, index)

h = 1;

dist_right_left= nan(1,ceil(size_y(2)-size_y(1)));
dist_top_bottom= nan(1,ceil(size_x(2)-size_x(1))); 

for pt_x = [size_x(1), size_x(2)]
    for pt_y = size_y(1):size_y(2)
        dist_right_left(1,h) = Distance(index(1),index(2), pt_x , pt_y); 
        h = h+1;
    end
end

h = 1;

for pt_y = [size_y(1), size_y(2)]
    for pt_x= size_x(1):size_x(2)
        dist_top_bottom(1,h)= Distance(index(1),index(2), pt_x , pt_y);
        h= h+1;
    end
end

all_dist = union(dist_top_bottom, dist_right_left);
dist = min(all_dist);

full_size_x= size_x(2)- size_x(1);
full_size_y= size_y(2)- size_y(1);

%divide by arena width to normalize
norm_top_bottom= dist_top_bottom/full_size_x;
norm_right_left= dist_right_left/full_size_y;

all_norm_dist = union(norm_top_bottom, norm_right_left);
norm_dist = min(all_norm_dist);

disp('');

