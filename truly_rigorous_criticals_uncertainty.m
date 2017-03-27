TR = min(T)/TC_fit;
% TR = 0.7

s = 1;

A_range = linspace(A_low,A_high,30);
b_range = linspace(b_low,b_high,30);
TC_range = linspace(TC_low,TC_high,40);
pc_range = linspace(pc_low,pc_high,40);

PDF = ones(length(A_range)*length(b_range)*length(pc_range)*length(TC_range),1);
TC = 0 * PDF;
rhoC = TC;
PC = TC;
ZC = TC;

for g=1:length(A_range)
    
    for h=1:length(b_range)

        for i=1:length(pc_range)
    
            for j=1:length(TC_range)
      
       SE = ((y - (pc_range(i) + A_range(g)*(TC_range(j)-T)))./sy).^2 + ((z - (b_range(h)*(TC_range(j)-T).^beta))./sz).^2;
 
       SSE = sum(SE);
       
       PDF(s) = fcdf((SSE-SSE_fit)/(sigma*p),p,n-p); % I assume that sigma = sigma2 in this code
       
       TC(s) = TC_range(j);
       
       rhoC(s) = pc_range(i);
       
       PC(s) = PC_from_Rackett_D(TC_range(j),pc_range(i),MW/pc_range(i),A_range(g),b_range(h),TR,D);
       
       ZC(s) = PC(s) * MW / (TC_range(j) * pc_range(i) * 8.314472);
       
       s=s+1;
       
            end
            
       end
       
   end
    
end

ii = 20000;

TC_norm = zeros(ii,1);
rhoC_norm = TC_norm;
PC_norm = TC_norm;
ZC_norm = TC_norm;

for i =1:ii
    
r=fcdf(frnd(p,n-p),p,n-p);

dif=abs(r-PDF);

[dif_min,I] = min(dif);

TC_norm(i) = TC(I);

rhoC_norm(i) = rhoC(I);

PC_norm(i) = PC(I);

ZC_norm(i) = ZC(I);

end

[TC_count, TC_centers] = hist(TC_norm,100);

figure
hist(TC_norm,100)

[rhoC_count, rhoC_centers] = hist(rhoC_norm,100);

figure
hist(rhoC_norm,100)

[PC_count, PC_centers] = hist(PC_norm,100);

figure
hist(PC_norm,100)

[ZC_count, ZC_centers] = hist(ZC_norm,100);

figure
hist(ZC_norm,100)

cd ../Parameter_space

[TC_min, TC_max, true_conf_TC] = integrate_histogram(TC_count,TC_centers,alpha);

[rhoC_min, rhoC_max, true_conf_rhoC] = integrate_histogram(rhoC_count,rhoC_centers,alpha);

[PC_min, PC_max, true_conf_PC] = integrate_histogram(PC_count,PC_centers,alpha);

[ZC_min, ZC_max, true_conf_ZC] = integrate_histogram(ZC_count,ZC_centers,alpha);

cd ../GEMC_Simulations_Analysis

% rhoC_min/pc_fit
% rhoC_max/pc_fit
% TC_min/TC_fit
% TC_max/TC_fit
% (PC_max-PC_min)/(PC_max+PC_min)
% (ZC_max-ZC_min)/(ZC_max+ZC_min)