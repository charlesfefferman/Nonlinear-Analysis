clear

% C2 from Martin's example file in Towhee

% T = [178 197 217 236 256 275]; 
% pv = [0.0021766 0.0056274 0.0099620 0.019879 0.030961 0.055561];
% pL = [0.55108 0.52686 0.49902 0.46898 0.43228 0.39554];
% sigmaL = [0.00087 0.0015355 0.0011414 0.0017548 0.0052414 0.0052225];
% sigmav = [0.00011795 0.00051137 0.00031008 0.0016763 0.0044995 0.0079753];
% 
% % Best fit parameters
% 
% A_fit = 5.6111576*10^-4;
% b_fit = 0.1168528;
% pc_fit = 0.2059202;
% TC_fit = 304.0091191;
% 
% % Ranges used
% 
% % Ranges used when sigmay=sigmaz
% % A_range = 0.00042:0.000005:0.0007;
% % b_range = 0.113:0.0001:0.1198;
% % pc_range = 0.19:0.0003:0.221;
% % TC_range = 296:0.2:316;
% 
% % Ranges used with the incorrect propagation of error
% A_range = 0.00048:0.0000025:0.00064;
% b_range = 0.1136:0.0001:0.1198;
% pc_range = 0.1957:0.0002:0.2153;
% TC_range = 296.4:0.2:313.4;

% C8_exact, all the 390,440,490,515 data with uncertainties

T = [390 390	390	390	390	390	390	390	390	390	390	390	390	390	390	390	390	390	390	390	390	390	440	440	440	440	440	440	440	440	440	440	440	440	440	440	440	440	440	440	440	440	440	440	490	490	490	490	490	490	490	490	490	490	490	490	490	490	490	490	490	490	490	490	490	490	515	515	515	515	515	515	515	515	515	515	515	515	515	515	515	515	515	515	515	515	515	515];
pv = [0.004327939	0.004184564	0.003917159	0.004101543	0.004477834	0.003507394	0.004370386	0.004295049	0.003909465	0.004095172	0.003961186	0.003628007	0.003696285	0.00372563	0.004538258	0.003791495	0.004025541	0.003892735	0.003682306	0.004455524	0.003812552	0.004219012	0.012571479	0.012331145	0.012285679	0.012126438	0.012875404	0.012683835	0.012679334	0.012977049	0.012062284	0.012652116	0.012510796	0.011913201	0.011516719	0.012130097	0.013270889	0.012597948	0.012811975	0.011882749	0.012281456	0.013343433	0.011970981	0.012817599	0.032072304	0.031102697	0.033176657	0.031307972	0.032473293	0.032543163	0.033180546	0.03079688	0.032405674	0.031062345	0.033426393	0.033677784	0.03288743	0.03282621	0.030675448	0.030836496	0.031768803	0.030750877	0.033692475	0.032863545	0.031560556	0.031015949	0.05051073	0.050830246	0.04899213	0.057110752	0.052651101	0.052339284	0.049552693	0.053307964	0.050235571	0.052363289	0.052616845	0.054453852	0.053871097	0.049989063	0.052184507	0.048926774	0.049999652	0.049765002	0.049955562	0.052363289	0.053323735	0.05073475];
pL = [0.62588879	0.624029464	0.62480371	0.624295007	0.626001224	0.623824273	0.624325807	0.625056439	0.624467903	0.625360839	0.62459378	0.624175523	0.624014819	0.624714229	0.626538872	0.624632971	0.624956424	0.625354585	0.624977206	0.625403293	0.624078045	0.625875774	0.573141775	0.573362883	0.573233316	0.572219544	0.574339626	0.574855147	0.575721706	0.57317294	0.573599977	0.572957209	0.573499894	0.572745119	0.573631613	0.572926795	0.573486909	0.574746432	0.574074923	0.57238545	0.573332909	0.573911377	0.572891604	0.574369315	0.5149841	0.511010445	0.512082774	0.512010323	0.511868857	0.515346328	0.511909917	0.511482209	0.513196845	0.511338411	0.513673407	0.51400598	0.513282354	0.512709949	0.511393191	0.511791603	0.511763086	0.513327369	0.512734432	0.511936205	0.512849102	0.5114893	0.475287289	0.474352627	0.470077928	0.476141671	0.47379976	0.473303127	0.475469186	0.475587343	0.477135568	0.47685196	0.475542761	0.474066296	0.477039424	0.476481347	0.476874839	0.47430863	0.475447446	0.469115182	0.47350821	0.47685196	0.476022927	0.475198579];
% sigmay = [0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000484676	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000567086	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.000946223	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146	0.001716146];
% sigmaz = [0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.000581324	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.00075617	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.001167651	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394	0.002302394];
% 
% % Best fit parameters
% 
A_fit = 4.191675*10^-4;
b_fit = 0.1177868;
pc_fit = 0.2385685;
TC_fit = 570.6878129;
% 
% % Ranges used
% 
A_range = 0.00041:0.000001:0.000428;
b_range = 0.1174:0.00002:0.1181;
pc_range = 0.2372:0.00005:0.2399;
TC_range = 569.75:0.05:571.6;

% C48 literature values from NERD

T = [780 810 840 860 870];
pv = [0.004 0.007 0.012 0.015 0.018];
pL = [0.538 0.490 0.442 0.404 0.382];
sigmaL = [0.006 0.006 0.005 0.012 0.009];
sigmav = [0.002 0.002 0.003 0.004 0.003];

% Best fit parameters

% A_fit = 7.6854753*10^-4;
% b_fit = 0.1120167;
% pc_fit = 0.1729306;
% TC_fit = 908.3512869;

% Best fit parameters when not using y and z, just pL and pv

% A_fit = 7.4462474*10^-4;
% b_fit = 0.1109362;
% pc_fit = 0.172611;
% TC_fit = 911.1614661;

% Best fit parameters when old y and z model was used

% A_fit = 7.6849734*10^-4;
% b_fit = 0.112367;
% pc_fit = 0.1728521;
% TC_fit = 908.0854856;

% Best fit parameters when new y and z model was used

A_fit = 7.6812344*10^-4;
b_fit = 0.1122975;
pc_fit = 0.1727803;
TC_fit = 908.2266955;

% Old, using just 6,8,10
% ey = 8.729809627769269;
% ez = 6.22;
% pey = 6.672793638522566*10^-7;
% pez = (3.45*10^-6)*2;

% New uncertainty model
ey = 9.052;
ez = 5.59;
pey = 5.4*10^-7;
pez = (5.92*10^-6)*2;

% Ranges used

% Ranges used when not propagation y and z, just using pL and pv
% A_range = -0.002:0.0001:0.002;
% b_range = 0:0.002:0.18;
% pc_range = 0:0.01:0.28;
% TC_range = 880:20:1800;

% Ranges used when old y and z model was used at the 95% confidence level
% A_range = -0.0001:0.000025:0.00165;
% b_range = 0.097:0.0005:0.129;
% % pc_range = 0.26:0.0003125:0.2625;
% % TC_range = 895.75:0.125:904.5;
% pc_range = 0.0825:0.00125:0.2625;
% TC_range = 882.25:0.25:952;

% Ranges used when new y and z model was used at the 95% confidence level
A_range = -0.0002:0.0001:0.00175;
b_range = 0.097:0.002:0.129;
% pc_range = 0.26:0.0003125:0.2625;
% TC_range = 895.75:0.125:904.5;
% pc_range = 0.0725:0.0025:0.2725;
% TC_range = 883.25:1:953.25;

% Using 2*data for n

A_range = (6.2:0.1:9.1)*10^-4;
b_range = 0.109:0.00025:0.115;
pc_range = 0.1575:0.00025:0.1875;
TC_range = 903:0.1:914;

% Ranges used when sigmay=sigmaz
% A_range = ;
% b_range = ;
% pc_range = ;
% TC_range = ;

% Ranges used with the incorrect propagation of error
% Refined range
% A_range = 0.0002:0.0000125:0.00135;
% b_range = 0.082:0.0005:0.13;
% pc_range = 0.055:0.00125:0.149;
% TC_range = 1000:0.25:1030;

% Full range
% A_range = 0.0002:0.000025:0.00135;
% b_range = 0.082:0.001:0.13;
% pc_range = 0.055:0.00125:0.23;
% TC_range = 880:0.25:1030;



n = 2*length(T);

beta = 0.32;

p = 4;

y = (pv + pL)/2;
z = pL - pv;

% sigmay = ((sigmaL.^2 + sigmav.^2).^(0.5))/2;
% sigmaz = (sigmaL.^2 + sigmav.^2).^(0.5);

% SE_fit = (y - (pc_fit + A_fit*(TC_fit-T))).^2 + (z - (b_fit*(TC_fit-T).^beta)).^2;

% SE_fit = ((y - (pc_fit + A_fit*(TC_fit-T)))./sigmay).^2 + ((z - (b_fit*(TC_fit-T).^beta))./sigmaz).^2;

% SE_fit = ((pL - (pc_fit + A_fit*(TC_fit-T) + 1/2*b_fit*(TC_fit-T).^beta))./sigmaL).^2 + ((pv - (pc_fit + A_fit*(TC_fit-T) - 1/2*b_fit*(TC_fit-T).^beta))./sigmav).^2;

SE_fit = ((y - (pc_fit + A_fit*(TC_fit-T)))./(pey*exp(ey*T/TC_fit))).^2 + ((z - (b_fit*(TC_fit-T).^beta))./(pez*exp(ez*T/TC_fit))).^2;

% Include the estimate of the uncertainties

% SE_fit = SE_fit./((sigmayz.^2));

SSE_fit = sum(SE_fit);

sigma = SSE_fit/(n-p);

RHS = sigma * (n + p * (finv(0.95^p,p,n-p)-1));

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
      
%        SE = (y - (pc_range(i) + A_range(g)*(TC_range(j)-T))).^2 + (z - (b_range(h)*(TC_range(j)-T).^beta)).^2;
 
%        SE = ((y - (pc_range(i) + A_range(g)*(TC_range(j)-T)))./sigmay).^2 + ((z - (b_range(h)*(TC_range(j)-T).^beta))./sigmaz).^2;
 
%        SE = ((pL - (pc_range(i) + A_range(g)*(TC_range(j)-T) + 1/2*b_range(h)*(TC_range(j)-T).^beta))./sigmaL).^2 + ((pv - (pc_range(i) + A_range(g)*(TC_range(j)-T) - 1/2*b_range(h)*(TC_range(j)-T).^beta))./sigmav).^2;
 
        SE = ((y - (pc_range(i) + A_range(g)*(TC_range(j)-T)))./(pey*exp(ey*T/TC_range(j)))).^2 + ((z - (b_range(h)*(TC_range(j)-T).^beta))./(pez*exp(ez*T/TC_range(j)))).^2;
 
%        SE = SE./((sigmayz.^2));

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

% subplot(3,2,1)
% plot(A_ext,b_ext)
% subplot(3,2,2)
% plot(pc_ext,b_ext)
% subplot(3,2,3)
% plot(A_ext,TC_ext)
% subplot(3,2,4)
% plot(pc_ext,TC_ext)
% subplot(3,2,5)
% plot(A_ext,pc_ext)
% subplot(3,2,6)
% plot(b_ext,TC_ext)

pc_high-pc_low
TC_high-TC_low

% hold
% 
% contour(TC_range,pc_range,confidence_region)
% 
% hold
% 
% % To plot just the pc and TC region

pc_space = 0.00125;

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

hold

plot(pc_ext,TC_ext)
plot(pc_plot,TC_upper);
plot(pc_plot,TC_lower);

hold