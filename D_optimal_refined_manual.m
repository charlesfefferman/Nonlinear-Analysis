clear pc_plot TC_lower TC_upper TC_ext pc_ext A_ext b_ext

beta = 0.32;

[TC_fit, pc_fit, A_fit, b_fit, sy, sz] = towhee_error_model(T,pv,pL);

n = 2*length(T);

p = 4;

y = (pv + pL)/2;
z = (pL - pv)/2;

SE_fit = ((y - (pc_fit + A_fit*(TC_fit-T)))./sy).^2 + ((z - (b_fit*(TC_fit-T).^beta))./sz).^2; % This includes uncertainties

SSE_fit = sum(SE_fit);

sigma = SSE_fit/(n-p);

RHS = sigma * (n + p * (finv(0.95^p,p,n-p)-1)); % No 0.95^p because I want the true 95%

% For C24 with D-optimal 200
% A_range = linspace(0.978,1.023,30)*A_fit;
% b_range = linspace(0.996,1.004,30)*b_fit;
% pc_range = linspace(0.9884,1.0116,400)*pc_fit;
% TC_range = linspace(0.9971,1.004,400)*TC_fit;

% For C24 with 600,680,715 and 1600 molecules
% A_range = linspace(0.976,1.024,30)*A_fit;
% b_range = linspace(0.9966,1.0034,30)*b_fit;
% pc_range = linspace(0.9906,1.0094,100)*pc_fit;
% TC_range = linspace(0.9979,1.0022,100)*TC_fit;

% For C48 for Nath et al.
% A_range = linspace(0.8102,1.1926,30)*A_fit;
% b_range = linspace(0.975,1.0217,30)*b_fit;
% pc_range = linspace(0.8920,1.1022,400)*pc_fit;
% TC_range = linspace(0.9919,1.0095,400)*TC_fit;

% % For C48 with Exp-6 700,800, as originally done in paper before 700
% finished
% A_range = linspace(0.8102,1.1926,30)*A_fit;
% b_range = linspace(0.975,1.0217,30)*b_fit;
% pc_range = linspace(0.8920,1.1022,50)*pc_fit;
% TC_range = linspace(0.9919,1.0095,50)*TC_fit;

% % For C48 with Exp-6 700,800
% A_range = linspace(0.9562,1.0438,30)*A_fit;
% b_range = linspace(0.9923,1.0077,30)*b_fit;
% pc_range = linspace(0.9781,1.0218,200)*pc_fit;
% TC_range = linspace(0.995,1.0052,200)*TC_fit;

% For C36 for TraPPE 650 775
% A_range = linspace(0.985,1.015,30)*A_fit;
% b_range = linspace(0.9975,1.0025,30)*b_fit;
% pc_range = linspace(0.992,1.008,200)*pc_fit;
% TC_range = linspace(0.998,1.002,200)*TC_fit;

% For C36 TraPPE 1600 730,775,810
% A_range = linspace(0.9609,1.0391,30)*A_fit;
% b_range = linspace(0.9955,1.0045,30)*b_fit;
% pc_range = linspace(0.9840,1.0144,200)*pc_fit;
% TC_range = linspace(0.9974,1.0026,200)*TC_fit;

% For C36 for Nath et al.
% A_range = linspace(0.8089,1.1911,30)*A_fit;
% b_range = linspace(0.9796,1.0191,30)*b_fit;
% pc_range = linspace(0.9098,1.0917,200)*pc_fit;
% TC_range = linspace(0.9927,1.009,200)*TC_fit;

% For C36 for NERD 650 765
% A_range = linspace(0.9778,1.0222,30)*A_fit;
% b_range = linspace(0.9965,1.0035,30)*b_fit;
% pc_range = linspace(0.9893,1.0107,200)*pc_fit;
% TC_range = linspace(0.997,1.003,200)*TC_fit;

% For C36 for NERD 1600 730 765 800 (old)
% A_range = linspace(0.9501,1.0501,30)*A_fit;
% b_range = linspace(0.9943,1.0057,30)*b_fit;
% pc_range = linspace(0.982,1.018,50)*pc_fit;
% TC_range = linspace(0.9965,1.0035,50)*TC_fit;

% For C36 for NERD 1600 730 765 800 (new)
% A_range = linspace(0.958,1.042,30)*A_fit;
% b_range = linspace(0.995,1.005,30)*b_fit;
% pc_range = linspace(0.9865,1.0135,200)*pc_fit;
% TC_range = linspace(0.9975,1.0025,200)*TC_fit;

% For C48 for Exp-6 with 700,800,835 and 1600 molecules
% A_range = linspace(0.94,1.06,30)*A_fit;
% b_range = linspace(0.99,1.01,30)*b_fit;
% pc_range = linspace(0.97,1.03,30)*pc_fit;
% TC_range = linspace(0.992,1.008,30)*TC_fit;

% For C48 for Exp-6 with 700,800 and 1600 molecules
% A_range = linspace(0.932,1.068,30)*A_fit;
% b_range = linspace(0.987,1.012,30)*b_fit;
% pc_range = linspace(0.9645,1.035,30)*pc_fit;
% TC_range = linspace(0.9928,1.008,30)*TC_fit;

% For C48 for NERD with 830, 860 and 1600 molecules
% A_range = linspace(0.79,1.21,40)*A_fit;
% b_range = linspace(0.98,1.02,30)*b_fit;
% pc_range = linspace(0.948,1.052,100)*pc_fit;
% TC_range = linspace(0.993,1.007,100)*TC_fit;

% For C48 for NERD with 800 (single run), 830, 860 and 1600 molecules
% A_range = linspace(0.885,1.115,40)*A_fit;
% b_range = linspace(0.9875,1.0125,30)*b_fit;
% pc_range = linspace(0.968,1.032,200)*pc_fit;
% TC_range = linspace(0.995,1.005,200)*TC_fit;

% For C48 for TraPPE with 730, 800, 830, 860 and 1600 molecules
% A_range = linspace(0.978,1.022,30)*A_fit;
% b_range = linspace(0.9965,1.0035,30)*b_fit;
% pc_range = linspace(0.99,1.01,30)*pc_fit;
% TC_range = linspace(0.998,1.002,30)*TC_fit;

% For C48 for TraPPE with 800, 830, 860 and 1600 molecules
A_range = linspace(0.925,1.075,40)*A_fit;
b_range = linspace(0.991,1.009,40)*b_fit;
pc_range = linspace(0.975,1.025,200)*pc_fit;
TC_range = linspace(0.995,1.005,200)*TC_fit;

s=1;
for g=1:length(A_range)
for h=1:length(b_range)
for i=1:length(pc_range)
for j=1:length(TC_range)
SE = ((y - (pc_range(i) + A_range(g)*(TC_range(j)-T)))./sy).^2 + ((z - (b_range(h)*(TC_range(j)-T).^beta))./sz).^2;
SSE = sum(SE);
if SSE < RHS
A_ext(s) = A_range(g);
b_ext(s) = b_range(h);
pc_ext(s) = pc_range(i);
TC_ext(s) = TC_range(j);
s=s+1;
if SSE < SSE_fit % This is to make sure that we actually have the global minimum. If not, rerun the analysis with the new optimum (after plugging in as a new Mathcad guess)
new_best_fit = [A_range(g) b_range(h) pc_range(i) TC_range(j)];
SSE_fit = SSE;
end
end
end
end
end
end
% This eliminates the superfluous zero elements
A_ext = A_ext(:,1:(s-1));
b_ext = b_ext(:,1:(s-1));
pc_ext = pc_ext(:,1:(s-1));
TC_ext = TC_ext(:,1:(s-1));
A_low_temp = min(A_ext);
A_high_temp = max(A_ext);
b_low_temp = min(b_ext);
b_high_temp = max(b_ext);
pc_low_temp = min(pc_ext);
pc_high_temp = max(pc_ext);
TC_low_temp = min(TC_ext);
TC_high_temp = max(TC_ext);

figure
subplot(3,2,1)
plot(A_ext,b_ext)
subplot(3,2,2)
plot(pc_ext,b_ext)
subplot(3,2,3)
plot(A_ext,TC_ext)
subplot(3,2,4)
plot(pc_ext,TC_ext)
subplot(3,2,5)
plot(A_ext,pc_ext)
subplot(3,2,6)
plot(b_ext,TC_ext)

% % To plot just the pc and TC region

pc_scan=pc_range;

k=1;

for h=1:length(pc_scan)
    
j=1;

    for i=1:length(pc_ext)
    
	if pc_ext(i) == pc_scan(h)

%         if (pc_ext(i) - pc_scan(h)) < 0.00000001
    
    TC_scan(j) = TC_ext(i);

    j=j+1;

        end

    end

% I did it this way because, for some reason, not every pc in the range of 
% pc_scan actually has an accepted point.  Somehow I missed the A,B and TC 
% combinations that were required.

    if j>1
       
    TC_upper(k) = max(TC_scan);
    TC_lower(k) = min(TC_scan);
    pc_plot(k) = pc_scan(h);
    
    k=k+1;
    
    end
    
TC_scan=TC_fit;

end

figure
hold

% scatter(pc_ext,TC_ext,'x')
plot(pc_plot,TC_upper);
plot(pc_plot,TC_lower);

hold

TC_low=min(TC_lower);

min(A_ext)/A_fit
max(A_ext)/A_fit
min(b_ext)/b_fit
max(b_ext)/b_fit
min(pc_ext)/pc_fit
max(pc_ext)/pc_fit
min(TC_ext)/TC_fit
max(TC_ext)/TC_fit
TC_high=max(TC_upper);
