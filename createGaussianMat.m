function [gaussian_mat]= createGaussianMat(mat, centers, centervalue)

gsize= size(mat);
new_size=gsize+20;
%kernalsize= sigma;

% for r=:kernalsize
%     for c=1:kernalsize
%         for m=1:length(centers)
%         gaussian_mat(r,c) = gaussC(r,c, sigma, centers(m,:));
%         end
%     end
% end

% if mod(parms.sigma,2)~= 0
%     parms.sigma= parms.sigma+1;
% end

val=centervalue;
%padded_pre_gaussian_mat= PreSmoothGaussian(gsize, centers, val);
%parms.sigma=2;
%gaussian_mat2= SmoothGaussian(padded_pre_gaussian_mat, parms);
%gaussian_mat2= gaussian_mat2(11:new_size(1)-10, 11:new_size(2)-10);


 %padded_pre_gaussian_mat= PreSmoothGaussian(gsize, centers, val*3.2);
 %parms.sigma=4;
 %gaussian_mat3= SmoothGaussian(padded_pre_gaussian_mat, parms);
 %gaussian_mat= gaussian_mat3(11:new_size(1)-10, 11:new_size(2)-10);
% 
% 
padded_pre_gaussian_mat= PreSmoothGaussian(gsize, centers, val*0.56);
 parms.sigma=1.5;
 gaussian_mat1= SmoothGaussian(padded_pre_gaussian_mat, parms);
 gaussian_mat= gaussian_mat1(11:new_size(1)-10, 11:new_size(2)-10);

disp('');

%plotPlot(mat, centervalue)


% function val = gaussC(x, y, sigma, center)
% xc = center(1);
% yc = center(2);
% exponent = ((x-xc).^2 + (y-yc).^2)./(2*sigma);
% val       = (exp(-exponent));
%  end

disp('');

function padded_pre_gaussian_mat = PreSmoothGaussian(gsize, centers, centervalue)

pre_gaussian_mat= zeros(gsize);

for h= 1:length(centers)
pre_gaussian_mat(centers(h,1), centers(h,2))= centervalue*25.25;
end

padded_pre_gaussian_mat= zeros(gsize + 20); % pad with edges
new_size= size(padded_pre_gaussian_mat);
padded_pre_gaussian_mat(11:new_size(1)-10, 11:new_size(2)-10)= ...
    pre_gaussian_mat;

function plotPlot(gm, mean_val)

imagesc(gm);axis image; axis off;
title(sprintf('%0.1f Hz', mean_val));
