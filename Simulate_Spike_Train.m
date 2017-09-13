function [spk_t]=Simulate_Spike_Train(pos_x,pos_y,pos_t,rate_map,parms,num)

dt=0.001;
max_t=max(pos_t);
max_x = ceil(max(pos_x));
max_y = ceil(max(pos_y));
min_x = min(floor(pos_x));
min_y = min(floor(pos_y));

% divid the environment into spatial bins
axis_x = min_x:parms.bin_size:max_x;
axis_y = min_y:parms.bin_size:max_y;

rnd=rand(num,length(0:dt:max_t));
rate=zeros(1,length(0:dt:max_t));
prob= zeros(1,length(0:dt:max_t));
count2=1;

spk_t= nan(num,1000); %preallocate for speed

size_rm= size(rate_map);

lambda = 1;
epsilon = 10^-4; 
pos_neg = [-1 1];

for t=max_t:(-1*dt):0
    
    [~,i]=min(abs(t-pos_t)); %Old way
    %[~,i] = histc(t, pos_t); %Newer way
        
    % finding the rate and probabilty
    [~,x_ind] =  min(abs(pos_x(i)-axis_x));
    [~,y_ind] =  min(abs(pos_y(i)-axis_y));

    x_ind(x_ind>size_rm(2))=size_rm(2);
    y_ind(y_ind>size_rm(1))=size_rm(1);
    
    rate(count2)=rate_map(y_ind,x_ind);
    prob(count2)=rate(count2)*dt;
    
    lambda = lambda * (1 + (epsilon * pos_neg(round(rand) + 1)));
    
    for h=1:num 
        if rnd(h,count2) < (prob(count2) * lambda) % && t-spk_t(end)>parms.refractoryPeriod
            
            ind=find(isnan(spk_t(h,:)));
            
            if isempty(ind) %if reached end of preallocated spk array
               spk_t(:,end+1)=nan; 
               ind=length(spk_t);
            end    
            
            ind=ind(1);
            spk_t(h,ind)=t;
        end
    end
    
    count2=count2+1;
    
end    

disp('')

