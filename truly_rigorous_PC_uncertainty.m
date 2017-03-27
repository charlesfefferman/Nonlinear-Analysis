TR = min(T)/TC_fit;
% TR = 0.7

s = 1;

% A_range = linspace(min(A_range),max(A_range),30);
% b_range = linspace(min(b_range),max(b_range),30);
% pc_range = linspace(min(pc_range),max(pc_range),40);
% TC_range = linspace(min(TC_range),max(TC_range),40);

PDF = ones(length(A_range)*length(b_range)*length(pc_range)*length(TC_range),1);
PC = 0 * PDF;
ZC = PC;

for g=1:length(A_range)
    
    for h=1:length(b_range)

        for i=1:length(pc_range)
    
            for j=1:length(TC_range)
      
       SE = ((y - (pc_range(i) + A_range(g)*(TC_range(j)-T)))./sy).^2 + ((z - (b_range(h)*(TC_range(j)-T).^beta))./sz).^2;
 
       SSE = sum(SE);
       
       PDF(s) = fcdf((SSE-SSE_fit)/(sigma*p),p,n-p); % I assume that sigma = sigma2 in this code
       
       PC(s) = PC_from_Rackett_D(TC_range(j),pc_range(i),MW/pc_range(i),A_range(g),b_range(h),TR,D);
       
       ZC(s) = PC(s) * MW / (TC_range(j) * pc_range(i) * 8.314472);
       
       s=s+1;
       
            end
            
       end
       
   end
    
end

for i =1:65500
    
r=fcdf(frnd(p,n-p),p,n-p);

dif=abs(r-PDF);

[dif_min,I] = min(dif);

PC_norm(i) = PC(I);

ZC_norm(i) = ZC(I);

end

[PC_count, PC_centers] = hist(PC_norm,1000);

figure
hist(PC_norm,50)

[ZC_count, ZC_centers] = hist(ZC_norm,1000);

figure
hist(ZC_norm,50)

cd ../Parameter_space

[PC_min, PC_max, true_conf_PC] = integrate_histogram(PC_count,PC_centers,alpha);

[ZC_min, ZC_max, true_conf_ZC] = integrate_histogram(ZC_count,ZC_centers,alpha);

cd ../GEMC_Simulations_Analysis

(PC_max-PC_min)/(PC_max+PC_min)
(ZC_max-ZC_min)/(ZC_max+ZC_min)