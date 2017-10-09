function out_mat=Cross_Correlation(mat1,mat2)

 [ma,na] = size(mat1);
 [mb,nb] = size(mat2);
 mc = max([ma+mb-1,ma,mb]);
 nc = max([na+nb-1,na,nb]);
 
 out_mat = nan(mc,nc);
 

i_size = size(mat2,1); j_size = size(mat2,2);
[work_mat,npad_i,npad_j] = pad_edges(mat1,mat2,1);

   for i = 1:size(out_mat,1)
        for j = 1:size(out_mat,2)
                
        % for each i and j, choose the correct sub-mat (size of mat 2) to
        % multiply with mat2
        
        sub_mat = work_mat(npad_i+i-floor(i_size):npad_i+i-1, ...
        npad_j+j-floor(j_size):npad_j+j-1  ); 
        nan_sub_mat=sub_mat .* mat2;                                             
        notnan_inds = find(~isnan(nan_sub_mat));  %normalized to the number of nontnan components (average)
        
        n=length(notnan_inds);
        
        if n < 20
            out_mat(i,j) = NaN;
            continue;
        end
        
        sigma_x_y =sum(nan_sub_mat(notnan_inds));
        sigma_x =      sum(sub_mat(notnan_inds));
        sigma_y =      sum(mat2(notnan_inds));
        sigma_x2 = sum(sub_mat(notnan_inds).^2);
        sigma_y2 = sum(mat2(notnan_inds).^2);
        
        out_mat(i,j) = (n*sigma_x_y - sigma_x.*sigma_y) ./ ...
                                    sqrt(n*sigma_x2-sigma_x.^2) ./ ...
                                    sqrt(n*sigma_y2-sigma_y.^2);
        

         end % for j
    end % for i
    
disp('')

function out_mat = nanconv2(mat,h)

out_mat = mat;
 nan_mat = isnan(mat);
 
 % dilate nan_mat
 
SE = strel('disk', 2);
nan_mat =  ~imdilate(~nan_mat,SE);

i_size = size(h,1); j_size = size(h,2);
[work_mat,npad_i,npad_j] = pad_edges(mat,h,2);
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

%figure; imagesc(out_mat);

disp('')
 
function [out_mat,npad_i,npad_j] = pad_edges(mat,h,l)

npad_ij = ceil(size(h)/l);
npad_i = npad_ij(1);
npad_j = npad_ij(2);
in_size = size(mat);
out_size = in_size + [2*npad_i 2*npad_j];
out_mat = nan(out_size);
out_mat(npad_i+1:npad_i+in_size(1),npad_j+1:npad_j+in_size(2)) = mat;

disp('')
