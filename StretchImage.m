function [stretched_image] = StretchImage(original_image, original_size, new_size)


factor_i= new_size(1)/original_size(1);
factor_j= new_size(2)/original_size(2);
tform = maketform('affine',[factor_j 0 0; 0 factor_i 0; 0 0 1]);
stretched_image = imtransform(original_image,tform);