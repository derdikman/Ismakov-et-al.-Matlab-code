function rate_mat=SmoothGaussian(rate_mat,parms) 
 %create window
 sigma=parms.sigma;
 size_h= [13 13]; 
  h=fspecial('gaussian',size_h,sigma);
  
  %figure;imagesc(h)
 %smooth the rate mat
rate_mat=nanconv2(rate_mat,h);
%figure;imagesc(rate_mat);title('rate map');axis equal
