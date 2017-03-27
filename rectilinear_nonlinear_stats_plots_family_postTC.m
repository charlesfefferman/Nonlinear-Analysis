clear

% tetracosane C24
 T1 = [580 590 600 610 625 675 700 725 750];
 pv1 = [0.0010309 0.0015445 0.0021023 0.0025044 0.003415 0.009588 0.016419 0.026868 0.032001];
 pL1 = [0.611756 0.6060484 0.5909036 0.5878981 0.57159 0.52602 0.50403 0.47456 0.43912];

% pentacosane C25
 T2 = [580 590 600 610 615 650 690 730 750];
 pv2 = [0.00109 0.001412 0.001406 0.002087415 0.00265 0.004574594 0.010672267 0.018401342 0.031339834];
 pL2 = [0.616996 0.608843 0.59835 0.592602248 0.587282 0.555748238 0.519263592 0.471308627 0.45664106];

% hexacosane C26
 T3 = [610 655 700 740 760];
 pv3 = [ 0.001423	0.00429041	0.010017451	0.023451055	0.026356036];
 pL3  = [ 0.592659	0.557560998	0.516726592	0.477174855	0.441255922];
% 615 0.001808302 0.59482208 data excluded b/c inconsistent with new data

% heptacosane C27
 T4 = [580 610 620 660 700 745 765];
 pv4 = [0.00036 0.001345 0.15333E-02 0.45482E-02  0.91630E-02 0.19332E-01 0.24939E-01];
 pL4 = [0.620732 0.602083 0.58842 0.55523 0.52438 0.47300 0.44285];
% The 580 data is questionable b/c not equil, but it improved regression

% octacosane C28
 T5 = [580 600 610 625 665 710 750 770];
 pv5 = [0.000304 0.000655 0.001123 0.17462E-02 0.48022E-02 0.11550E-01 0.17752E-01 0.25336E-01];
 pL5 = [0.624137 0.612879 0.602442 0.59089 0.56211 0.51428 0.46550 0.45160];

 
A1 = 4.01008 * 10^-4;
b1 = 0.106109;
pc1 = 0.2117755;
TC1 = 816.1717;

A2 = 4.00359 * 10^-4;
b2 = 0.105801;
pc2 = 0.2106519;
TC2 = 824.9692;

A3 = 4.01609 * 10^-4;
b3 = 0.105353;
pc3 = 0.2087054;
TC3 = 833.4005;

A4 = 4.07275 * 10^-4;
b4 = 0.104729;
pc4 = 0.205853;
TC4 = 841.4931;

A5 = 4.04992 * 10^-4;
b5 = 0.1043;
pc5 = 0.2053387;
TC5 = 849.2707;

n1 = length(T1);
n2 = length(T2);
n3 = length(T3);
n4 = length(T4);
n5 = length(T5);

beta = 0.32;
p = 3; % Only 3 now because we have fixed TC so that it satisfies the linear equation perfectly.

y1 = (pv1 + pL1)/2;
z1 = pL1 - pv1;

y2 = (pv2 + pL2)/2;
z2 = pL2 - pv2;

y3 = (pv3 + pL3)/2;
z3 = pL3 - pv3;

y4 = (pv4 + pL4)/2;
z4 = pL4 - pv4;

y5 = (pv5 + pL5)/2;
z5 = pL5 - pv5;

SE_fit1 = (y1 - (pc1 + A1*(TC1-T1))).^2 + (z1 - (b1*(TC1-T1).^beta)).^2;
SE_fit2 = (y2 - (pc2 + A2*(TC2-T2))).^2 + (z2 - (b2*(TC2-T2).^beta)).^2;
SE_fit3 = (y3 - (pc3 + A3*(TC3-T3))).^2 + (z3 - (b3*(TC3-T3).^beta)).^2;
SE_fit4 = (y4 - (pc4 + A4*(TC4-T4))).^2 + (z4 - (b4*(TC4-T4).^beta)).^2;
SE_fit5 = (y5 - (pc5 + A5*(TC5-T5))).^2 + (z5 - (b5*(TC5-T5).^beta)).^2;

SSE_fit1 = sum(SE_fit1);
SSE_fit2 = sum(SE_fit2);
SSE_fit3 = sum(SE_fit3);
SSE_fit4 = sum(SE_fit4);
SSE_fit5 = sum(SE_fit5);

sigma1 = SSE_fit1/(n1-p);
sigma2 = SSE_fit2/(n2-p);
sigma3 = SSE_fit3/(n3-p);
sigma4 = SSE_fit4/(n4-p);
sigma5 = SSE_fit5/(n5-p);

RHS1 = sigma1 * (n1 + p * (finv(0.95^p,p,n1-p)-1));
RHS2 = sigma2 * (n2 + p * (finv(0.95^p,p,n2-p)-1));
RHS3 = sigma3 * (n3 + p * (finv(0.95^p,p,n3-p)-1));
RHS4 = sigma4 * (n4 + p * (finv(0.95^p,p,n4-p)-1));
RHS5 = sigma5 * (n5 + p * (finv(0.95^p,p,n5-p)-1));

% These are supposed to be slightly larger than the extrema found below, just so that you can verify this range.
% Be careful, for the b_range I realized that the displayed number is
% rounded off, so make sure it really found the min or max

% Ranges used for 24
% A_range = 0.000338:0.000002:0.000468; 
% b_range = 0.10432:0.00004:0.10774;
% pc_range = 0.2:0.0002:0.225;
% TC_range = 808:0.2:825;

% Ranges used for 25
% A_range = 0.000346:0.000002:0.000454; 
% b_range = 0.10392:0.00004:0.10694;
% pc_range = 0.197:0.0003:0.221;
% TC_range = 820:0.2:835;

% Ranges used for 26 (these data are worthless) Essentially you must have
% more data points to lock down pc-TC at all

% A_range = -0.0002:0.00002:0.001; 
% b_range = 0.08:0.0005:0.12;
% pc_range = 0.08:0.004:0.31;
% TC_range = 780:2:960;

% Ranges used for 27

%  A_range = 0.0003:0.000004:0.00052; 
%  b_range = 0.101:0.0001:0.108;
%  pc_range = 0.18:0.0005:0.23;
%  TC_range = 825:0.5:860;
 
% Ranges used for 28

A_range = 0.00031:0.000004:0.0005; 
b_range = 0.101:0.0001:0.11;
pc_range = 0.184:0.0001:0.226;

% No longer a range for TC

s=1;
t=1;
u=1;

for g=1:length(A_range)
    
    for h=1:length(b_range)

        for i=1:length(pc_range)
                   
       SE = (y1 - (pc_range(i) + A_range(g)*(TC1-T1))).^2 + (z1 - (b_range(h)*(TC1-T1).^beta)).^2;

       SSE = sum(SE);
       
       if SSE > RHS1 
           
       else
           
           %We start confidence_region at one value and change it only if
           %the SSE is less than the RHS.  This is because we want to plot
           %any pc-TC that is acceptable. Otherwise we would be
           %overwritting accepted points whenever one is rejected.
                
           A_ext(s) = A_range(g);
           b_ext(t) = b_range(h);
           pc_ext(u) = pc_range(i);
           
           s=s+1;
           t=t+1;
           u=u+1;
       end
       
            
       end
       
   end
    
end

A_low = min(A_ext);
A_high = max(A_ext);
b_low = min(b_ext);
b_high = max(b_ext);
pc_low = min(pc_ext)
pc_high = max(pc_ext)