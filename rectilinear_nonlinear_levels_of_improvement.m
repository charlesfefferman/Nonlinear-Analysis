clear

% This is the TraPPE ethane data
% I am trying to find an example where the CI becomes non-physical
% If I restrict TraPPE's data to just the three higher temperatures it
% overlaps at the 95% confidence level using the traditional approach. I
% want to see how it looks with the rigorous approach.

T = [178 197 217 236 256 275];
pv = [0.0021766 0.0056274 0.0099620 0.019879 0.030961 0.055561];
pL = [0.55108 0.52686 0.49902 0.46898 0.43228 0.39554];
sv = [0.00011795 0.00051137 0.00031008 0.0016763 0.0044995 0.0079753];
sL = [0.00087617 0.0015355 0.0011414 0.0017548 0.0052414 0.0052225];

m = 4;
n = length(T);

T = T(m:n);
pv = pv(m:n);
pL = pL(m:n);
sv = sv(m:n);
sL = sL(m:n);

sigmay = ((sL.^2 + sv.^2).^(0.5))/2;
sigmaz = sigmay;

% A_fit = 5.3211073*10^-4;
% b_fit = 0.058425;
% pc_fit = 0.2085556;
% TC_fit = 303.187101;

% This with m=4 and alpha =0.95
% A_range = (2.9:0.1:7.8)*10^-4;
% b_range = (5.5:0.02:6.1)*10^-2;
% pc_range = 0.191:0.0005:0.225;
% TC_range = 295.5:0.25:314.5;

% This with m=4 and alpha =0.995 and TT error
% A_range = (-2.6:0.2:13.2)*10^-4;
% b_range = (4.5:0.1:6.6)*10^-2;
% pc_range = 0.132:0.002:0.26;
% TC_range = 285:1:375;

% This with m=4 and alpha =0.995 and RR error
% A_range = (-3.8:0.4:15.4)*10^-4;
% b_range = (5.05:0.05:6.4)*10^-2;
% pc_range = 0.134:0.002:0.268;
% TC_range = 289:0.5:337.5;
% 
% % This is the same error model that RRT used (rigorous-rigorous-traditional)
% pea = 8.12195682246548*10^-9;
% pes = 1.9346962419217932*10^-8;
% ea = 13.880645592227626;
% es = 12.18250766453814;
% 
% sigmay = pea*exp(ea*T/TC_fit);
% sigmaz = pes*exp(es*T/TC_fit);

% This is with constant term in error model
A_fit = 5.6689432*10^-4;
b_fit = 0.0584767;
pc_fit = 0.2061589;
TC_fit = 303.0566329;

% This with m=4 and alpha =0.995 and RR error with constant term
% A_range = (-3.8:0.5:15.2)*10^-4;
% b_range = (5:0.05:6.5)*10^-2;
% pc_range = 0.142:0.001:0.262;
% TC_range = 289:0.25:335;

% With 0.995^2
% A_range = linspace(-8,18,40)*10^-4;
% b_range = linspace(4,7,40)*10^-2;
% pc_range = linspace(0.1,0.3,40);
% TC_range = linspace(280,380,40);

% With n=5*n0 and 0.995^2

% A_range = linspace(5,6.5,40)*10^-4;
% b_range = linspace(5.5,6,40)*10^-2;
% pc_range = linspace(0.195,0.215,40);
% TC_range = linspace(298,305,40);

% With n=5*3 and 0.995^2 and NERD's error

A_range = linspace(4,7,40)*10^-4;
b_range = linspace(5,6,40)*10^-2;
pc_range = linspace(0.18,0.22,50);
TC_range = linspace(290,320,50);
% 
% % This is the same error model that RRT used but with a constant term included (rigorous-rigorous-traditional)
b0a = 5*10^-4;
b1a = 2.25*10^-12;
b2a = 22.475;
b0s = 3.9875*10^-4;
b1s = 2.225*10^-14;
b2s = 26.55;

sigmay = b0a + b1a*exp(b2a*T/TC_fit);
sigmaz = b0s + b1s*exp(b2s*T/TC_fit);

% C8 data

% T = [390 440 490 515 543];
% pv = [0.004219 0.013343 0.032864 0.05073475 0.0840558];
% pL = [0.625876 0.573911 0.511936 0.475198579 0.418986443];
% % 
% % % Traditional length, traditional propagation of error:
% sv = [0.000625 0.002313 0.005560 0.008574752 0.013940418]; % These are for a single run, the traditional way
% sL = [0.003440 0.003881 0.004561 0.009607261 0.018147998];
% 
% sigmay = ((sL.^2 + sv.^2).^(0.5))/2;
% sigmaz = (sL.^2 + sv.^2).^(0.5)/2;
% 
% A_fit = 4.2279208*10^-4;
% b_fit = 0.0590918;
% pc_fit = 0.2391457;
% TC_fit = 569.3150885;
% 
% A_range = 0.000412:0.000002:0.000434;
% b_range = 0.0588:0.00004:0.0594;
% pc_range = 0.2375:0.00005:0.241;
% TC_range = 568:0.05:570.7;
% 
% % With 0.95^1
% 
% A_range = linspace(0.0004, 0.00044, 60);
% b_range = linspace(0.0585, 0.0595, 60);
% pc_range = linspace(0.236, 0.242, 60);
% TC_range = linspace(566,573,60);

% Traditional length, rigorous propagation of error:
% sa = [0.001007704 0.0024603 0.0043773 0.0086883 0.013896621];	% Found using traditional length of simulation
% ss = [0.001047761 0.0020383 0.0025878 0.0027253 0.008290142];
% 
% sigmay = sa;
% sigmaz = ss;
% 
% A_fit = 4.2412298*10^-4;
% b_fit = 0.0590773;
% pc_fit = 0.2389539;
% TC_fit = 569.3593042;
% 
% A_range = 0.000414:0.000001:0.000434;
% b_range = 0.0589:0.00002:0.0593;
% pc_range = 0.2373:0.00005:0.2406;
% TC_range = 568.6:0.025:570.2;
% 
% % With 0.95^1
% 
% A_range = linspace(0.0004, 0.00044, 60);
% b_range = linspace(0.0585, 0.0595, 60);
% pc_range = linspace(0.236, 0.242, 60);
% TC_range = linspace(568,571,60);


% Rigorous length, traditional propagation of error:
% sigmav = [0.000300663 0.000463719 0.001039459 0.002035368 0.005824739]; % These are the extended standard deviations for liquid and vapor for C8
% sigmaL = [0.000740536 0.000845006 0.00117979 0.002115684 0.00559524];
% 
% sigmay = ((sigmaL.^2 + sigmav.^2).^(0.5))/2;
% sigmaz = (sigmaL.^2 + sigmav.^2).^(0.5)/2;
% 
% A_fit = 4.2361004*10^-4;
% b_fit = 0.0590759;
% pc_fit = 0.2389505;
% TC_fit = 569.4509494;
% % 
% A_range = 0.000414:0.0000005:0.000433;
% b_range = 0.0588:0.00002:0.0594;
% pc_range = 0.2375:0.00005:0.2404;
% TC_range = 568.1:0.05:570.9;
% 
% % With 0.95^1
% 
% A_range = linspace(0.0004, 0.00044, 60);
% b_range = linspace(0.0585, 0.0595, 60);
% pc_range = linspace(0.236, 0.242, 60);
% TC_range = linspace(566,573,60);

% Rigorous length, rigorous propagation of error:
% sigmaa = [0.000484676 0.000567086 0.000946223 0.001727084 0.004908074];
% sigmas = [0.000290662 0.000378085 0.000583825 0.001151781 0.002628723];
% 
% sigmay = sigmaa;
% sigmaz = sigmas;
% 
% A_fit = 4.2353772*10^-4;
% b_fit = 0.0590783;
% pc_fit = 0.2389733;
% TC_fit = 569.4180815;
% 
% A_range = 0.000412:0.0000005:0.000435;
% b_range = 0.0588:0.00001:0.0593;
% pc_range = 0.2373:0.00005:0.2407;
% TC_range = 568.4:0.025:570.5;
% 
% % With 0.95^1
% 
% A_range = linspace(0.0004, 0.00044, 60);
% b_range = linspace(0.0585, 0.0595, 60);
% pc_range = linspace(0.236, 0.242, 60);
% TC_range = linspace(568,571,60);

n = 5*length(T);

beta = 0.32;

p = 4;

alpha = 0.995;

y = (pv + pL)/2;
z = (pL - pv)/2;

SE_fit = ((y - (pc_fit + A_fit*(TC_fit-T)))./sigmay).^2 + ((z - (b_fit*(TC_fit-T).^beta))./sigmaz).^2;

SSE_fit = sum(SE_fit);

sigma = SSE_fit/(n-p);

RHS = sigma * (n + p * (finv(alpha^2,p,n-p)-1));

% confidence_region = ones(length(pc_range),length(TC_range));

Ext = zeros(1,(length(A_range)*length(b_range)*length(pc_range)*length(TC_range)));
A_ext = Ext;
b_ext = Ext;
pc_ext = Ext;
TC_ext = Ext;
clear Ext;

s=1;

for g=1:length(A_range)
    
    for h=1:length(b_range)

        for i=1:length(pc_range)
    
            for j=1:length(TC_range)
      
 
       SE = ((y - (pc_range(i) + A_range(g)*(TC_range(j)-T)))./sigmay).^2 + ((z - (b_range(h)*(TC_range(j)-T).^beta))./sigmaz).^2;
 
       SSE = sum(SE);
       
       if SSE < RHS 
                     
           %We start confidence_region at one value and change it only if
           %the SSE is less than the RHS.  This is because we want to plot
           %any pc-TC that is acceptable. Otherwise we would be
           %overwritting accepted points whenever one is rejected.
           
%            confidence_region(i,j) = 0; % We only want to plot pc vs TC
           
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

% AbpcTC = [A_ext; b_ext; pc_ext; TC_ext]';

A_low = min(A_ext);
A_high = max(A_ext);
b_low = min(b_ext);
b_high = max(b_ext);
pc_low = min(pc_ext);
pc_high = max(pc_ext);
TC_low = min(TC_ext);
TC_high = max(TC_ext);

figure
subplot(3,2,1)
plot(A_ext,b_ext)
subplot(3,2,2)
plot(pc_ext,b_ext)
subplot(3,2,3)
plot(A_ext,TC_ext)
subplot(3,2,4)
hold
plot(pc_ext,TC_ext)
scatter(pc_fit,TC_fit,'r')
hold
subplot(3,2,5)
plot(A_ext,pc_ext)
subplot(3,2,6)
plot(b_ext,TC_ext)

(pc_high-pc_low)/2
(TC_high-TC_low)/2

% hold
% 
% contour(TC_range,pc_range,confidence_region)
% 
% hold
% 
% % To plot just the pc and TC region

pc_space = 0.0001;

pc_scan=min(pc_range):pc_space:max(pc_range);

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

% figure
% hold
% 
% plot(pc_ext,TC_ext)
% plot(pc_plot,TC_upper);
% plot(pc_plot,TC_lower);
% 
% hold

A_low
A_high
b_low
b_high
pc_low
pc_high
TC_low
TC_high