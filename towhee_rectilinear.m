% clear

% This is the NERD ethane literature data, testing out towhee
%Towhee gave: 
% rhoc = 0.199325 0.00125
% TC = 310.8197 11.4413

% T = [175.4 213.8 252.2 271 280.8];
% rhog = [0.005 0.009 0.023 0.041 0.049];
% rhol = [0.543 0.496 0.438 0.405 0.381];
% errg = [0.001 0.001 0.002 0.002 0.003];
% errl = [0.004 0.003 0.004 0.005 0.007];

% This is the TraPPE ethane data provided in towhee
% Towhee gave:
% rhoc 0.2059 0.00624
% TC 303.98 5.692

T = [178 197 217 236 256 275];
rhog = [0.0021766 0.0056274 0.0099620 0.019879 0.030961 0.055561];
rhol = [0.55108 0.52686 0.49902 0.46898 0.43228 0.39554];
errg = [0.00011795 0.00051137 0.00031008 0.0016763 0.0044995 0.0079753];
errl = [0.00087617 0.0015355 0.0011414 0.0017548 0.0052414 0.0052225];

% I am trying to find an example where the CI becomes non-physical
% If I restrict TraPPE's data to just the three higher temperatures it
% overlaps at the 95% confidence level.

m = 4;
n = length(T);

T = T(m:n);
rhog = rhog(m:n);
rhol = rhol(m:n);
errg = errg(m:n);
errl = errl(m:n);

% This is the TraPPE data for C6 that I simulated

% T = [460 480];
% rhog = [0.062596 0.112918609];
%     rhol = [0.464863 0.425490519];
%     errg = [0.015925734 0.03639983];
%     errl = [0.009046722 0.016039826];

% This is the NERD 48 data

% T = [780 810 840 860 870];
% rhog = [0.004 0.007 0.012 0.015 0.018];
% rhol = [0.538 0.490 0.442 0.404 0.382];
% errl = [0.006 0.006 0.005 0.012 0.009];
% errg = [0.002 0.002 0.003 0.004 0.003];

n = 2*length(T);

beta = 0.32;

rhoa = (rhog + rhol) / 2;
erra = sqrt(errg.^2 + errl.^2)/2;

delrho = (rhol - rhog) / 2;
rhos = delrho;
errs = erra;
rhos = rhos.^(1/beta);
errs = rhos * (1/beta) .* errs ./ delrho;

[as, bs, sigas, sigbs] = towhee_fit(T,rhos,errs);

TC = -as/bs;
dTC = TC * sqrt((sigas/as)^2 + (sigbs/bs)^2);
a0 = as^beta;
da = ((as^beta)*beta*sigas/as);

[ar, br, sigar, sigbr] = towhee_fit(T,rhoa,erra);

b1 = -br*TC;
db = b1*sqrt((sigbr/br)^2+(dTC/TC)^2);
rhoc = ar-b1;
drhoc = sqrt(sigar^2+db^2);

if (n-4) == 4
    
    tstat = 2.776445;
    
elseif (n-4) == 6
    
    tstat = 2.446912;
        
elseif (n-4) == 8
    
    tstat = 2.306004;
    
else
    
    tstat = 1.3*tinv(0.95,n-4);
    
end

% tstat = 14.09; % This is for the graph of traditional vs rigorous

CI = tstat/sqrt(n); % The 1.3 is to approximately make it the two-tailed solution in the range of 4-6 temperatures

TC_low = TC-CI*dTC;
TC_high = TC+CI*dTC;
rhoc_low = rhoc-CI*drhoc;
rhoc_high = rhoc+CI*drhoc;