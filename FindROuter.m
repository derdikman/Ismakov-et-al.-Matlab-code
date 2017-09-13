function R_outer=FindROuter(acorr)
    
  % calculate all the extrema points in the spatial autocorrelation 
[zmax,imax,zmin,imin]= Extrema2(acorr);
[i,j]=ind2sub(size(acorr),imax);


%put all extrema points in dist
dist(:,1)=j;
dist(:,2)=i;
n=length(i);

%calculate the distance of all extrema to the central peak and put them in
%column 3
dist(1:n,3)=sqrt(  (i(1:n)-i(1)).^2 + (j(1:n)-j(1)).^2);
 
 % sort the hexonal peaks by distance to the centeral peak
 [score,ind]=sort(dist(:,3));
 dist=dist(ind,:);
 %zmax=zmax(ind);
 R=dist(2,3);
  count=1;
 i=2;
 hex_peaks(1,:,:)=dist(1,:,:);
 
% finds the first 6 closest peaks to the central peak
while count<7 && i<=length(dist)

    % calculate the min distance of point i from all other already chosen
    % points
    min_dist_peaks=min(sqrt(  (hex_peaks(1:size(hex_peaks,1),1)-dist(i,1)).^2 +...
     (hex_peaks(1:size(hex_peaks,1),2)-dist(i,2)) .^2));

    % point i needs to be on the right side (we choos only half the point cause its semetrical)
    % and the distance of point i from all other already chosen points
    % needs to be higher than R/2
    if dist(i,1)>=dist(1,1) && min_dist_peaks>(R/1.5)
        hex_peaks(count,:,:)=dist(i,:,:); 
        count=count+1; 
        hex_peaks(count,1)=dist(1,1)-dist(i,1)+dist(1,1);
        hex_peaks(count,2)=dist(1,2)-dist(i,2)+dist(1,2);

        count=count+1;
    end    
    i=i+1;
end
 
R_outer=max(hex_peaks(:,3))*1.15;  
    
