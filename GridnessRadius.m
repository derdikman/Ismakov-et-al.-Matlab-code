function [gridness2] = GridnessRadius...
        (org_Acor,R_outer,image_id)

    
   %% compute distance from center
[val,center_y]=max(max(org_Acor));
[val,center_x]=max(max(org_Acor'));
[Y,X] = ndgrid(1:1:size(org_Acor,1), 1:1:size(org_Acor,2));
dist_from_center=sqrt((Y-center_y).^2+(X-center_x).^2);

%fig_image = figure;imagesc(dist_from_center);

%cd(parms.dir_save_data);
%saveas(fig_image,sprintf('distance_from_center_i=%d.fig',image_id));
%saveas(fig_image,sprintf('distance_from_center_i=%d.jpg',image_id));


%% making sure that outer radius of the ring (R_outer) is not bigger than the distance matrix
R_outer=min([min(dist_from_center(1,:)),min(dist_from_center(:,1)),...
    min(dist_from_center(size(dist_from_center,1),:))...
min(dist_from_center(:,size(dist_from_center,2))),R_outer]);

%% compute inner radius of the anulus (ring)
R_inner=ceil(min(dist_from_center(org_Acor<0.1)));

    
%extract the original anulus (ring) from Acor
org_Ring=org_Acor(dist_from_center<=R_outer & dist_from_center>=R_inner);


%% plot original ring
%  figure;
%  tmp=org_Acor;
% tmp(dist_from_center>R_outer | dist_from_center<R_inner)=nan
%  imagesc(tmp)
 

% make sure that after rotation and interpulation the center will remain the maximum point.
    org_Acor(center_x,center_y)=10;

    for jj = 2:6
         
       %% rotate the auto-correlation 
       rot_Acor=imrotate(org_Acor,(jj-1)*30,'bicubic');
                              
         %% compute distance from new center
        [val,tmp_center_x]=max(max(rot_Acor));
        [val,tmp_center_y]=max(max(rot_Acor'));
        [Y,X] = ndgrid(1:1:size(rot_Acor,1), 1:1:size(rot_Acor,2));
        tmp_dist_from_center=sqrt((Y-tmp_center_y).^2+(X-tmp_center_x).^2);
        
        % extract the anulus(ring)
        rot_Ring=rot_Acor(tmp_dist_from_center<=R_outer & tmp_dist_from_center>=R_inner);
        
%         figure;
%          tmp_rot1=rot_Acor
%          tmp_rot1(tmp_dist_from_center>R_outer | tmp_dist_from_center<R_inner)=nan
%           imagesc(tmp_rot1)
%          
        if length(rot_Ring)~=length(org_Ring) 
            gridness2=nan;
            return 
        end
        
        %% compute pearson correlation between rotate Acor and original Acor
        corrValues(jj) = PointCorr(org_Ring,rot_Ring);
        
       clear rot_Ring tmp_center_x tmp_center_y tmp_dist_from_center Y X
        
    end
    
     

%    % min of higher correlation at 60 and 120 degree rotation
     min_rot_60_120 = min([corrValues(3),corrValues(5)]);
    % max of lower correlations 30,90,150 rotation
    max_rot_30_90_150 = max([corrValues(2),corrValues(4),corrValues(6)]);
    
    %% calculate gridness min(60,120)-max(30,90,150)
    gridness2 = min_rot_60_120 - max_rot_30_90_150;

%% different way to calculate gridness
     gridness1=mean(([corrValues(3),corrValues(5)]))-...
         mean([corrValues(2),corrValues(4),corrValues(6)]);
