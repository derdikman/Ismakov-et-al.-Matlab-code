function out_mat = nanconv2(mat,h)

out_mat = mat;
 nan_mat = isnan(mat);
 
 % dilate nan_mat
 
SE = strel('disk', 2);
nan_mat =  ~imdilate(~nan_mat,SE);

i_size = size(h,1); j_size = size(h,2);
[work_mat,npad_i,npad_j] = PadEdges(mat,h,2);
for i = 1:size(mat,1)
    for j = 1:size(mat,2)
        
        % for each i and j, choose the correct sub-mat (size of h) to multiply with h
        
        sub_mat = work_mat(npad_i+i-floor(i_size/2):npad_i+i+floor(i_size/2), ...
                                                          npad_j+j-floor(j_size/2):npad_j+j+floor(j_size/2)  ); % assumes h is odd in number
                                                      
        notnan_inds = find(~isnan(sub_mat));  
        
        if ~isempty(notnan_inds)
            sum_h = sum(h(notnan_inds));   % normalize to the places without a NaN
            out_mat(i,j) = nansum(nansum(sub_mat .* h));
            out_mat(i,j) = out_mat(i,j)/sum_h;
        end
        
    end % for j
end % for i

out_mat(nan_mat) = NaN;

disp('')
 