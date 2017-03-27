clear

% NERD literature

% C2

NC = 2;

T = [175.4 213.8 252.2 271 280.8];
pv = [0.005 0.009 0.023 0.041 0.049];
pL = [0.543 0.496 0.438 0.405 0.381];
sv = [0.001 0.001 0.002 0.002 0.003];
sL = [0.004 0.003 0.004 0.005 0.007];

% C3

NC = 3;

T = [250 281 312 340 344];
pv = [0.005 0.012 0.027 0.055 0.058];
pL = [0.546 0.516 0.474 0.428 0.418];
sv = [0.001 0.002 0.003 0.003 0.004];
sL = [0.006 0.005 0.007 0.007 0.01];

% % C4
% 
NC = 4;

T = [295 327 360 380 392];
pv = [0.005 0.014 0.024 0.040 0.051];
pL = [0.566 0.533 0.487 0.460 0.438];
sv = [0.001 0.002 0.004 0.003 0.005];
sL = [0.004 0.006 0.008 0.007 0.01];
% 
% % C5
% 
NC = 5;

T = [336 372 402 439];
pv = [0.006 0.013 0.023 0.055];
pL = [0.568 0.528 0.487 0.433];
sv = [0.002 0.002 0.003 0.004];
sL = [0.003 0.004 0.006 0.006];
% 
% % C6
% 
NC = 6;

T = [370 400 430 450 470];
pv = [0.0065 0.013 0.023 0.037 0.051];
pL = [0.577 0.546 0.506 0.480 0.439];
sv = [0.0005 0.001 0.003 0.003 0.006];
sL = [0.006 0.008 0.005 0.009 0.011];
% 
% % C8
% 
NC = 8;

T = [400 440 490 540];
pv = [0.003 0.009 0.026 0.063];
pL = [0.61 0.567 0.508 0.424];
sv = [0.001 0.001 0.004 0.004];
sL = [0.003 0.006 0.005 0.008];
% 
% % C10
% 
NC = 10;

T = [450 480 520 565 580];
pv = [0.005 0.011 0.023 0.050 0.060];
pL = [0.602 0.570 0.525 0.469 0.439];
sv = [0.002 0.002 0.003 0.006 0.007];
sL = [0.005 0.006 0.005 0.005 0.007];
% 
% % C12
% 
NC = 12;

T = [525 550 575 600 625];
pv = [0.009 0.013 0.023 0.034 0.053];
pL = [0.558 0.532 0.5 0.461 0.415];
sv = [0.001 0.002 0.003 0.004 0.005];
sL = [0.004 0.004 0.006 0.009 0.012];

% % C16 remember this isn't the best data b/c I had to steal it from graph
% 
NC = 16;

T = [520 600 650 680 690];
pv = [0.001 0.012 0.038 0.05 0.055];
pL = [0.61 0.53 0.45 0.42 0.41];
sv = [0.001 0.002 0.002 0.003 0.005]; % Guesses
sL = [0.004 0.004 0.006 0.006 0.008]; 

% C24

NC = 24;

T = [651.3 701.5 746.6 770];
pv = [0.006 0.017 0.032 0.048];
pL = [0.544 0.492 0.428 0.385];
sv = [0.002 0.002 0.004 0.005];
sL = [0.005 0.005 0.007 0.009];

% C36

% NC = 36;
% 
% T = [735 760 790 820];
% pv = [0.007 0.010 0.017 0.024];
% pL = [0.520 0.484 0.440 0.378];
% sv = [0.001 0.002 0.003 0.005];
% sL = [0.007 0.006 0.009 0.009];
% 
% % % C48
% % 
% NC = 48;
% 
% T = [780 810 840 860 870];
% pv = [0.004 0.007 0.012 0.015 0.018];
% pL = [0.538 0.490 0.442 0.404 0.382];
% sL = [0.006 0.006 0.005 0.012 0.009];
% sv = [0.002 0.002 0.003 0.004 0.003];

[TC, rhoc, dTC, drhoc] = towhee_regression(T,pv,pL,sv,sL);

n = length(T);
N = 20;

if n == 4
    
    tstat = 1.99085;
    
elseif n == 5
    
    tstat = 1.984467455;
           
else
    
    tstat = 1.3*tinv(0.95,N*n-2);
    
end

CI = tstat/sqrt(N); % The 1.3 is to approximately make it the two-tailed solution in the range of 4-6 temperatures

TC_low = TC-CI*dTC;
TC_high = TC+CI*dTC;
rhoc_low = rhoc-CI*drhoc;
rhoc_high = rhoc+CI*drhoc;

MW = 12.0107*NC+1.00794*(2*NC+2);

VC_high = MW/rhoc_low;
VC = MW/rhoc;
VC_low = MW/rhoc_high;

% subplot(2,1,1)
% hold
% scatter(NC,VC)
% plot(NC*[1,1],[VC_low,VC_high],'--')
% hold
% 
% subplot(2,1,2)
% hold
% scatter(NC,TC)
% plot(NC*[1,1],[TC_low,TC_high],'--')
% hold

[NC VC_low VC VC_high TC_low TC TC_high];

plot_data = [VC-VC_low, VC_high-VC, CI*dTC];

[NC rhoc drhoc/sqrt(N) dTC/sqrt(N)]