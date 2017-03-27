clear
clc

% % Validation data
% 
T = [178 197 217 256 275 279 283 288];
pv = [0.0023 0.0053 0.0111 0.0350 0.0598 0.0648 0.0739 0.090];
pL = [0.5512 0.5262 0.4984 0.4342 0.3937 0.3835 0.3726 0.3589];
sv = [0.0001 0.0001 0.0001 0.0003 0.0003 0.0005 0.0009 0.002];
sL = [0.0001 0.0001 0.0001 0.0002 0.0003 0.0004 0.0003 0.0006];

T = T(5:end);
pv = pv(5:end);
pL = pL(5:end);
sv = sv(5:end);
sL = sL(5:end);
T = [T T T T T T T T];
pv = [pv pv pv pv pv pv pv pv];
pL = [pL pL pL pL pL pL pL pL];
sv = [sv sv sv sv sv sv sv sv];
sL = [sL sL sL sL sL sL sL sL];

sy = sqrt(sv.^2 + sL.^2)/2;
sz = sy;

[TC_fit, pc_fit, dTC, dpc, A_fit, b_fit] = towhee_regression(T,pv,pL,sv,sL); %Make sure beta is set to 0.326

n = 2*length(T);

beta = 0.326;
p = 4;

y = (pv + pL)/2;
z = (pL - pv)/2;

SE_fit = ((y - (pc_fit + A_fit*(TC_fit-T)))./sy).^2 + ((z - (b_fit*(TC_fit-T).^beta))./sz).^2; % This includes uncertainties

SSE_fit = sum(SE_fit);

sigma = SSE_fit/(n-p);

RHS = sigma * (n + p * (finv(0.95^2,p,n-p)-1)); % To compare TC and rhoc
% For the growth part we start with a pretty coarse grid.  Once we have
% grown as much as we want, then we go to a finer grid.

A_range = linspace(0.788,1.22,100)*A_fit;
b_range = linspace(0.988,1.011,100)*b_fit;
pc_range = linspace(0.989,1.01,100)*pc_fit;
TC_range = linspace(0.997,1.004,100)*TC_fit;

% This is for the figures in my paper where two parameters are at the best
% fit values
A_range=A_fit;
b_range=b_fit;
% pc_range=pc_fit;
% TC_range=TC_fit;

% A_range = linspace(0.97,1.03,1000)*A_fit;
% b_range = linspace(0.998,1.002,1000)*b_fit;

% A_range = linspace(0.7797,1.221,1000)*A_fit;
% pc_range = linspace(0.9898,1.0102,1000)*pc_fit;

% A_range = linspace(0.9712,1.0287,1000)*A_fit;
% TC_range = linspace(0.9995,1.0005,1000)*TC_fit;

% b_range = linspace(0.99825,1.00175,1000)*b_fit;
% pc_range = linspace(0.9987,1.0013,1000)*pc_fit;

% b_range = linspace(0.9926,1.0074,1000)*b_fit;
% TC_range = linspace(0.9982,1.0019,1000)*TC_fit;
% 
pc_range = linspace(0.99865,1.0014,1000)*pc_fit;
TC_range = linspace(0.9995,1.0005,1000)*TC_fit;

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
      
        SE = ((y - (pc_range(i) + A_range(g)*(TC_range(j)-T)))./sy).^2 + ((z - (b_range(h)*(TC_range(j)-T).^beta))./sz).^2;
 
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

A_low/A_fit
A_high/A_fit
b_low/b_fit
b_high/b_fit
pc_low/pc_fit
pc_high/pc_fit
TC_low/TC_fit
TC_high/TC_fit

% subplot(3,2,1)
% plot(A_ext,b_ext)
% subplot(3,2,2)
% plot(b_ext,pc_ext)
% subplot(3,2,3)
% plot(A_ext,TC_ext)
% subplot(3,2,4)
plot(pc_ext,TC_ext)
% subplot(3,2,5)
% plot(A_ext,pc_ext)
% subplot(3,2,6)
% plot(b_ext,TC_ext)

% Ab_spacing = 10; % For speed, only 10 values will be used for A and b
% pcTC_spacing = 20; % pc and TC are more important
% iteration = 0;
% 
% % Converged is for the finer grid portion, growing is for the coarse grid
% converged = false;
% convergence = [false, false, false, false, false, false, false, false];
% growing = true;
% growth = [true,true,true,true,true,true,true,true];
% 
% % Since we are going to grow our confidence region we start with a very
% % small region.  These will grow with the coarse grid and then contract
% % back with the finer grid.
% 
% A_low = 0.99*A_fit;
% A_high = 1.01*A_fit;
% b_low = 0.995*b_fit;
% b_high = 1.005*b_fit;
% pc_low = 0.995*pc_fit;
% pc_high = 1.005*pc_fit;
% % pc_low = 0.2405;
% % pc_high = 0.24445;
% TC_low = 0.995*TC_fit;
% TC_high = 1.005*TC_fit;
% % TC_low = 616.75;
% % TC_high = 617.735;
% 
% 
% while converged ==false
%     
% A_spacing = (A_high-A_low)/(Ab_spacing-1);
% b_spacing = (b_high-b_low)/(Ab_spacing-1);
% pc_spacing = (pc_high-pc_low)/(pcTC_spacing-1);
% TC_spacing = (TC_high-TC_low)/(pcTC_spacing-1);
% 
% A_range = A_low:A_spacing:A_high;
% b_range = b_low:b_spacing:b_high;
% pc_range = pc_low:pc_spacing:pc_high;
% TC_range = TC_low:TC_spacing:TC_high;
% 
% Ext = zeros(1,(length(A_range)*length(b_range)*length(pc_range)*length(TC_range)));
% A_ext = Ext;
% b_ext = Ext;
% pc_ext = Ext;
% TC_ext = Ext;
% clear Ext
% 
% s=1;
% 
% for g=1:length(A_range)
%     
%     for h=1:length(b_range)
% 
%         for i=1:length(pc_range)
%     
%             for j=1:length(TC_range)
%       
%        SE = ((y - (pc_range(i) + A_range(g)*(TC_range(j)-T)))./sy).^2 + ((z - (b_range(h)*(TC_range(j)-T).^beta))./sz).^2;
%  
%        SSE = sum(SE);
%        
%        if SSE < RHS 
%                  
%            A_ext(s) = A_range(g);
%            b_ext(s) = b_range(h);
%            pc_ext(s) = pc_range(i);
%            TC_ext(s) = TC_range(j);
%            
%            s=s+1;
%            
%            if SSE < SSE_fit % This is to make sure that we actually have the global minimum. If not, rerun the analysis with the new optimum (after plugging in as a new Mathcad guess)
%                
%               new_best_fit = [A_range(g) b_range(h) pc_range(i) TC_range(j)]; 
%               SSE_fit = SSE;
%               
%            end
%            
%        end
%        
%             end
%             
%        end
%        
%    end
%     
% end
% 
% % This eliminates the superfluous zero elements
% A_ext = A_ext(:,1:(s-1));
% b_ext = b_ext(:,1:(s-1));
% pc_ext = pc_ext(:,1:(s-1));
% TC_ext = TC_ext(:,1:(s-1));
% 
% A_low_temp = min(A_ext);
% A_high_temp = max(A_ext);
% b_low_temp = min(b_ext);
% b_high_temp = max(b_ext);
% pc_low_temp = min(pc_ext);
% pc_high_temp = max(pc_ext);
% TC_low_temp = min(TC_ext);
% TC_high_temp = max(TC_ext);
% 
% if growing % The first part is to grow the regions
%     
%     if isempty(A_low_temp)
%         
%         factor = 10^(-(6+iteration));
%         low_factor = 1-factor;
%         high_factor = 1+factor;
%         
%         A_low = low_factor*A_fit;
%         A_high = high_factor*A_fit;
%         b_low = low_factor*b_fit;
%         b_high = high_factor*b_fit;
%         pc_low = low_factor*pc_fit;
%         pc_high = high_factor*pc_fit;
%         TC_low = low_factor*TC_fit;
%         TC_high = high_factor*TC_fit;
%         
%         
%     else
%     
%     if A_low_temp == A_low  
%     
%     A_low = A_low - 2*A_spacing;
%     growth(1) = true;
%     
%     else
%     
%     growth(1) = false;
%     
%     end
% 
%     if A_high_temp == A_high  
%     
%     A_high = A_high + 2*A_spacing;
%     growth(2) = true;
%       
%     else
%     
%     growth(2) = false;
%     
%     end
%     
%     if b_low_temp == b_low  
%     
%     b_low = b_low - 2*b_spacing;
%     growth(3) = true;
%     
%     else
%     
%     growth(3) = false;
%     
%     end
% 
%     if b_high_temp == b_high  
%     
%     b_high = b_high + 2*b_spacing;
%     growth(4) = true;
%     
%     else
%     
%     growth(4) = false;
%     
%     end
% 
%     if pc_low_temp == pc_low  
%     
%     pc_low = pc_low - 2*pc_spacing;
%     growth(5) = true;
%     
%     else
%     
%     growth(5) = false;
%     
%     end
% 
%     if pc_high_temp == pc_high  
%     
%     pc_high = pc_high + 2*pc_spacing;
%     growth(6) = true;
%     
%     else
%     
%     growth(6) = false;
%     
%     end
% 
%     if TC_low_temp == TC_low  
%     
%     TC_low = TC_low - 2*TC_spacing;
%     growth(7) = true;
%     
%     else
%     
%     growth(7) = false;
%     
%     end
% 
%     if TC_high_temp == TC_high  
%     
%     TC_high = TC_high + 2*TC_spacing;
%     growth(8) = true;
%     
%     else
%     
%     growth(8) = false;
%     
%     end
%     
%     end
%     
%     if growth(:)==false % Once all the regions incapsiluate the extrema we contract back down with a more refined grid
%         
%     growing=false;
%     Ab_spacing = 20; % For speed, only 20 values will be used for A and b
%     pcTC_spacing = 40; % pc and TC are more important
%     
%     end
%     
% else
% 
% if isempty(A_low_temp) % If one is empty they are all empty (could also use if s==1)
%     
%     A_low = A_low + A_spacing/2;
%     A_high = A_high - A_spacing/2;
%     b_low = b_low + b_spacing/2;
%     b_high = b_high - b_spacing/2;
%     pc_low = pc_low + pc_spacing/2;
%     pc_high = pc_high - pc_spacing/2;
%     TC_low = TC_low + TC_spacing/2;
%     TC_high = TC_high - TC_spacing/2;
%     
% else
% 
%     if A_low_temp == A_low  
%     
%    convergence(1) = true;  
%     
%     else
%     
%    A_low = A_low_temp - A_spacing/2;
%     
%     end
% 
%     if A_high_temp == A_high  
%     
%    convergence(2) = true;
%       
%     else
%     
%    A_high = A_high_temp + A_spacing/2;
%     
%     end
%     
%     if b_low_temp == b_low  
%     
%    convergence(3) = true;
%     
%     else
%     
%     b_low = b_low_temp - A_spacing/2;
%     
%     end
% 
%     if b_high_temp == b_high  
%     
%     convergence(4) = true;
%     
%     else
%     
%     b_high = b_high_temp + A_spacing/2;
%     
%     end
% 
%     if pc_low_temp == pc_low  
%     
%     convergence(5) = true;
%     
%     else
%     
%     pc_low = pc_low_temp - pc_spacing/2;
%     
%     end
% 
%     if pc_high_temp == pc_high  
%     
%     convergence(6) = true;
%     
%     else
%     
%     pc_high = pc_high_temp + pc_spacing/2;
%     
%     end
% 
%     if TC_low_temp == TC_low  
%     
%     convergence(7) = true;
%     
%     else
%     
%     TC_low = TC_low_temp - TC_spacing/2;
%     
%     end
% 
%     if TC_high_temp == TC_high  
%     
%     convergence(8) = true;
%     
%     else
%     
%     TC_high = TC_high_temp + TC_spacing/2;
%     
%     end
% 
% end
% 
% end
% 
% if convergence(:)==true % If all have converged we are done
%         
%     converged=true;
%     
% end
%     
% iteration = iteration + 1;
% 
% if iteration == 30 % To avoid infinite loops
%     
%     converged = true;
%     
% end
% 
% end
% 
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

% To plot just the pc and TC region

k=1;

for h=1:length(pc_range)
    
j=1;

    for i=1:length(pc_ext)
    
	if pc_ext(i) == pc_range(h)

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
    pc_plot(k) = pc_range(h);
    
    k=k+1;
    
    end
    
TC_scan=TC_fit;

end

figure
hold

% plot(pc_ext,TC_ext)
plot(pc_plot,TC_upper)
plot(pc_plot,TC_lower)

hold
% 
% TC_low=min(TC_lower);
% TC_high=max(TC_upper);