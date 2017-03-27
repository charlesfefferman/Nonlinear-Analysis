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


A = 4.031634e-4;

b1 = 0.1060823; 
b2 = 0.1054392; 
b3 = 0.1048676;
b4 = 0.1046103;
b5 = 0.1044448;

pc1 = 0.21137;
pc2 = 0.2094131;
pc3 = 0.2077101;
pc4 = 0.2063328;
pc5 = 0.2060522;

TC1 = 816.295217;
TC2 = 826.7979716;
TC3 = 835.3305909;
TC4 = 842.0696458;
TC5 = 848.270683;

n1 = length(T1);
n2 = length(T2);
n3 = length(T3);
n4 = length(T4);
n5 = length(T5);
n = n1+n2+n3+n4+n5;

beta = 0.32;
p = 16; % 16 now because all the parameters are regressed simultaneously

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

SE_fit1 = (y1 - (pc1 + A*(TC1-T1))).^2 + (z1 - (b1*(TC1-T1).^beta)).^2;
SE_fit2 = (y2 - (pc2 + A*(TC2-T2))).^2 + (z2 - (b2*(TC2-T2).^beta)).^2;
SE_fit3 = (y3 - (pc3 + A*(TC3-T3))).^2 + (z3 - (b3*(TC3-T3).^beta)).^2;
SE_fit4 = (y4 - (pc4 + A*(TC4-T4))).^2 + (z4 - (b4*(TC4-T4).^beta)).^2;
SE_fit5 = (y5 - (pc5 + A*(TC5-T5))).^2 + (z5 - (b5*(TC5-T5).^beta)).^2;

SSE_fit1 = sum(SE_fit1);
SSE_fit2 = sum(SE_fit2);
SSE_fit3 = sum(SE_fit3);
SSE_fit4 = sum(SE_fit4);
SSE_fit5 = sum(SE_fit5);

SSE_fit = SSE_fit1 + SSE_fit2 + SSE_fit3 + SSE_fit4 + SSE_fit5;

sigma2 = SSE_fit/(n-p);

RHS = sigma2 * (n + p * (finv(0.95^p,p,n-p)-1)); % Just one RHS because all one regression

% These are supposed to be slightly larger than the extrema found below, just so that you can verify this range.
% Be careful, for the b_range I realized that the displayed number is
% rounded off, so make sure it really found the min or max

A_range = 0.00038566:0.000002:0.00041966;
% b_range = 0; %0.0004 Failed approach
% pc_range = 0; %0.003
% TC_range = 0; %2

% Ranges used for 24
%  b_range1 = (b1-b_range):b_range:(b1+b_range);
%  pc_range1 = (pc1-pc_range):pc_range:(pc1+pc_range);
%  TC_range1 = (TC1-TC_range):TC_range:(TC1+TC_range);

% b_range1 = 0.1036:0.00008:0.1084;
% pc_range1 = 0.203:0.00025:0.219;
% TC_range1 = 806:0.25:828;

b_range5 = b5;
pc_range5=pc5;
TC_range5=TC5;
b_range1 = b1;
pc_range1=pc1;
TC_range1=TC1;
b_range3 = b3;
pc_range3=pc3;
TC_range3=TC3;
b_range4 = b4;
pc_range4=pc4;
TC_range4=TC4;

% Ranges used for 25
% b_range2 = (b2-b_range):b_range:(b2+b_range);
% pc_range2 = (pc2-pc_range):pc_range:(pc2+pc_range);
% TC_range2 = (TC2-TC_range):TC_range:(TC2+TC_range);

% Found by using the best fit values for all other parameters.
b_range2 = 0.103:0.0001:0.1078;
pc_range2 = 0.201:0.00025:0.217;
TC_range2 = 815:0.3:840;

% Ranges used for 26 (these data are worthless) Essentially you must have
% more data points to lock down pc-TC at all

% b_range3 = (b3-b_range):b_range:(b3+b_range);
% pc_range3 = (pc3-pc_range):pc_range:(pc3+pc_range);
% TC_range3 = (TC3-TC_range):TC_range:(TC3+TC_range);

% b_range3 = 0.1008:0.0001:0.1087;
% pc_range3 = 0.197:0.00025:0.217;
% TC_range3 = 820:0.5:855;

% Ranges used for 27

% b_range4 = (b4-b_range):b_range:(b4+b_range);
% pc_range4 = (pc4-pc_range):pc_range:(pc4+pc_range);
% TC_range4 = (TC4-TC_range):TC_range:(TC4+TC_range);

% b_range4 = 0.1018:0.0001:0.1073;
% pc_range4 = 0.197:0.0002:0.215;
% TC_range4 = 830:0.3:856;

% Ranges used for 28

% b_range5 = (b5-b_range):b_range:(b5+b_range);
% pc_range5 = (pc5-pc_range):pc_range:(pc5+pc_range);
% TC_range5 = (TC5-TC_range):TC_range:(TC5+TC_range);

% b_range5 = 0.10185:0.0001:0.1069;
% pc_range5 = 0.197:0.00025:0.214;
% TC_range5 = 837:0.3:863;

% confidence_region1 = ones(length(pc_range1),length(TC_range1));
% confidence_region2 = ones(length(pc_range2),length(TC_range2));
% confidence_region3 = ones(length(pc_range3),length(TC_range3));
% confidence_region4 = ones(length(pc_range4),length(TC_range4));
% confidence_region5 = ones(length(pc_range5),length(TC_range5));
Ext = zeros(1,(length(A_range)*length(b_range2)*length(pc_range2)*length(TC_range2)));
A_ext = Ext;
b_ext = Ext;
pc_ext = Ext;
TC_ext = Ext;
s=1;


for g=1:length(A_range) % Only one A range
    
    for h1=1:length(b_range1) 

        for i1=1:length(pc_range1)
    
            for j1=1:length(TC_range1)
                
                for h2=1:length(b_range2)
                    
                    for i2=1:length(pc_range2)
                        
                        for j2=1:length(TC_range2)
                            
                            for h3=1:length(b_range3)
                                
                                for i3=1:length(pc_range3)
                                    
                                    for j3=1:length(TC_range3)
                                        
                                        for h4=1:length(b_range4)
                                            
                                            for i4=1:length(pc_range4)
                                                
                                                for j4=1:length(TC_range4)
                                                    
                                                    for h5=1:length(b_range5)
                                                        
                                                        for i5=1:length(pc_range5)
                                                            
                                                            for j5=1:length(TC_range5)
      
       SE1 = (y1 - (pc_range1(i1) + A_range(g)*(TC_range1(j1)-T1))).^2 + (z1 - (b_range1(h1)*(TC_range1(j1)-T1).^beta)).^2;
       SE2 = (y2 - (pc_range2(i2) + A_range(g)*(TC_range2(j2)-T2))).^2 + (z2 - (b_range2(h2)*(TC_range2(j2)-T2).^beta)).^2;
       SE3 = (y3 - (pc_range3(i3) + A_range(g)*(TC_range3(j3)-T3))).^2 + (z3 - (b_range3(h3)*(TC_range3(j3)-T3).^beta)).^2;
       SE4 = (y4 - (pc_range4(i4) + A_range(g)*(TC_range4(j4)-T4))).^2 + (z4 - (b_range4(h4)*(TC_range4(j4)-T4).^beta)).^2;
       SE5 = (y5 - (pc_range5(i5) + A_range(g)*(TC_range5(j5)-T5))).^2 + (z5 - (b_range5(h5)*(TC_range5(j5)-T5).^beta)).^2;

       SSE1 = sum(SE1);
       SSE2 = sum(SE2);
       SSE3 = sum(SE3);
       SSE4 = sum(SE4);
       SSE5 = sum(SE5);

       SSE = SSE1 + SSE2 + SSE3 + SSE4 + SSE5;
       
       if SSE > RHS 
           
       else
           
           %We start confidence_region at one value and change it only if
           %the SSE is less than the RHS.  This is because we want to plot
           %any pc-TC that is acceptable. Otherwise we would be
           %overwritting accepted points whenever one is rejected.
           
%            confidence_region1(i1,j1) = 0; % We only want to plot pc vs TC
%            confidence_region2(i2,j2) = 0;
%            confidence_region3(i3,j3) = 0;
%            confidence_region4(i4,j4) = 0;
%            confidence_region5(i5,j5) = 0;
           
           A_ext(s) = A_range(g);
           b_ext(s) = b_range2(h2);
           pc_ext(s) = pc_range2(i2);
           TC_ext(s) = TC_range2(j2);
           s=s+1;
           
       end
                                                            end
                                                            
                                                        end
                                                        
                                                    end
                                                    
                                                end
                                                
                                            end
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
           end
            
       end
       
   end
    
end

% This eliminates the superfluous zero elements
A_ext(:,all(~A_ext,1))=[];
b_ext(:,all(~b_ext,1))=[];
pc_ext(:,all(~pc_ext,1))=[];
TC_ext(:,all(~TC_ext,1))=[];

AbpcTC = [A_ext; b_ext; pc_ext; TC_ext]';

A_low = min(A_ext);
A_high = max(A_ext);
b_low = min(b_ext);
b_high = max(b_ext);
pc_low = min(pc_ext);
pc_high = max(pc_ext);
TC_low = min(TC_ext);
TC_high = max(TC_ext);

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

% hold

% contour(TC_range1,pc_range1,confidence_region1)
% contour(TC_range2,pc_range2,confidence_region2)
% contour(TC_range3,pc_range3,confidence_region3)
% contour(TC_range4,pc_range4,confidence_region4)
% contour(TC_range5,pc_range5,confidence_region5)

% hold