clear

NC = [24 25 26 27 28];

TC = [816.171727 824.969164 833.400485 841.493145 849.270742];

y_fit = 1./TC;

x = 1 * NC.^(-0.5) + 1./(2*NC);

x_bar = mean(x);

y_bar_fit = mean(y_fit);

S_xx = sum((x - x_bar).^2);

S_xy_fit = sum((x - x_bar).*(y_fit - y_bar_fit));

m_fit = S_xy_fit / S_xx;

b_fit = y_bar_fit - m_fit * x_bar;

SSE_fit = sum((y_fit - (b_fit + m_fit * x)).^2 * 10^12);  % To avoid round-off errors

n = length(x);

p = 2;

sigma2 = SSE_fit / (n - p);

RHS = sigma2 * (n + p * (finv(0.95^p,p,n-p)-1));

% TC1_range = 809.5:1:824;
% TC2_range = 820.2:1:834;
% TC3_range = 800:5:870;
% TC4_range = 829:3:863;
% TC5_range = 836:2:858;

TC1_range = 816.16:0.0005:816.18;
TC2_range = 824.962:0.0005:824.974;
TC3_range = 833.3:0.01:833.5;
TC4_range = 841.489:0.0005:841.498;
TC5_range = 849.265:0.0005:849.28;

TC_region = zeros(1,5);
plot_region = zeros(length(TC1_range),length(TC2_range));
plot_region2 = zeros(length(TC4_range),length(TC5_range));

for f=1:length(TC1_range)
    
    for g=1:length(TC2_range)
        
        for h=1:length(TC3_range)
            
            for i=1:length(TC4_range)
                
                for j=1:length(TC5_range)
                    
                    TC_set = [TC1_range(f) TC2_range(g) TC3_range(h) TC4_range(i) TC5_range(j)];
                    
                    y = 1./TC_set;

                    y_bar = mean(y);

                    S_xy = sum((x - x_bar).*(y - y_bar));

                    m = S_xy / S_xx;

                    b = y_bar - m * x_bar;

                    SSE = sum((y - (b + m * x)).^2 * 10^12);  % To avoid round-off errors
 
                    
                    if SSE < RHS
                       
                        TC_region = [TC_region; TC_set];
                        plot_region(f,g) = 1;
                        plot_region2(i,j) = 1;
                                             
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

TC_region = TC_region(2:end,:);

subplot(2,1,1)
contour(TC2_range,TC1_range,plot_region)
subplot(2,1,2)
contour(TC5_range,TC4_range,plot_region2)