s = 1;

A_range = linspace(A_low,A_high,30);
b_range = linspace(b_low,b_high,30);
TC_range = linspace(TC_low,TC_high,40);
pc_range = linspace(pc_low,pc_high,40);

PDF = ones(length(A_range)*length(b_range)*length(pc_range)*length(TC_range),1);
TC = 0 * PDF;
rhoC = TC;

for g=1:length(A_range)
    
    for h=1:length(b_range)

        for i=1:length(pc_range)
    
            for j=1:length(TC_range)
      
       SE = ((y - (pc_range(i) + A_range(g)*(TC_range(j)-T)))./sy).^2 + ((z - (b_range(h)*(TC_range(j)-T).^beta))./sz).^2;
 
       SSE = sum(SE);
       
       PDF(s) = fcdf((SSE-SSE_fit)/(sigma*p),p,n-p); % I assume that sigma = sigma2 in this code
       
       TC(s) = TC_range(j);
       
       rhoC(s) = pc_range(i);
       
       s=s+1;
       
            end
            
       end
       
   end
    
end

ii = 10000;

TC_norm = zeros(ii,1);
rhoC_norm = TC_norm;

for i =1:ii
    
r=fcdf(frnd(p,n-p),p,n-p);

dif=abs(r-PDF);

[dif_min,I] = min(dif);

TC_norm(i) = TC(I);

rhoC_norm(i) = rhoC(I);

end

[TC_count, TC_centers] = hist(TC_norm,80);

figure
hist(TC_norm,80)

[rhoC_count, rhoC_centers] = hist(rhoC_norm,80);

figure
hist(rhoC_norm,80)

cd ../Parameter_space

[TC_min, TC_max, true_conf_TC] = integrate_histogram(TC_count,TC_centers,alpha);

[rhoC_min, rhoC_max, true_conf_rhoC] = integrate_histogram(rhoC_count,rhoC_centers,alpha);

cd ../GEMC_Simulations_Analysis

rhoC_min/pc_fit
rhoC_max/pc_fit
TC_min/TC_fit
TC_max/TC_fit