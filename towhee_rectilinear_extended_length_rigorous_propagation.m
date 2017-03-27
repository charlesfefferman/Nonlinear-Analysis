% clear

% Since this code uses the extended length we have used a model for the
% literature results.  However, for my actual simulation results I am not
% limited to this, if I want to I can find the standard deviations specific
% to C8, but the model might be more reliable.

% This is the NERD ethane literature data, testing out towhee
%Towhee gave: 
% rhoc = 0.199325 0.00125
% TC = 310.8197 11.4413

% TC = 310.8197;
% 
% T = [175.4 213.8 252.2 271 280.8];
% rhog = [0.005 0.009 0.023 0.041 0.049];
% rhol = [0.543 0.496 0.438 0.405 0.381];
% errg = [0.001 0.001 0.002 0.002 0.003];
% errl = [0.004 0.003 0.004 0.005 0.007];

% This is the TraPPE ethane data provided in towhee
% Towhee gave:
% rhoc 0.2059 0.00624
% TC 303.98 5.692

% TC = 303.98;
% 
T = [178 197 217 236 256 275];
rhog = [0.0021766 0.0056274 0.0099620 0.019879 0.030961 0.055561];
rhol = [0.55108 0.52686 0.49902 0.46898 0.43228 0.39554];
errg = [0.00011795 0.00051137 0.00031008 0.0016763 0.0044995 0.0079753];
errl = [0.00087617 0.0015355 0.0011414 0.0017548 0.0052414 0.0052225];
% 
% % I am trying to find an example where the CI becomes non-physical
% % If I restrict TraPPE's data to just the three higher temperatures it
% % overlaps at the 95% confidence level.
% 
m = 4;
n = length(T);

T = T(m:n);
rhog = rhog(m:n);
rhol = rhol(m:n);
errg = errg(m:n);
errl = errl(m:n);

n = 2*length(T);

% This is my TraPPE octane data
% Towhee gave:
% rhoc 0.239147 0.02396
% TC 569.312433 24.019290

% TC = 569;
% 
% T = [390 440 490 515 543];
% rhog = [0.004219 0.013343 0.032864 0.05073475 0.0840558];
% rhol = [0.625876 0.573911 0.511936 0.475198579 0.418986443];
% errg = [0.000300663 0.000463719 0.001039459 0.002035368 0.005824739]; % These are the extended standard deviations for liquid and vapor for C8
% errl = [0.000740536 0.000845006 0.00117979 0.002115684 0.00559524];
% erra = [0.000484676 0.000567086 0.000946223 0.001727084 0.004908074];
% errs = [0.000290662 0.000378085 0.000583825 0.001151781 0.002628723];
% 
beta = 0.32;

rhoa = (rhog + rhol) / 2;

delrho = (rhol - rhog) / 2;
rhos = delrho;
rhos = rhos.^(1/beta);

% If using the model, uncomment here and below
% Non-constant term error model
% pea = 8.12195682246548*10^-9;
% pes = 1.9346962419217932*10^-8;
% ea = 13.880645592227626;
% es = 12.18250766453814;

% Constant term error model
b0a = 5*10^-4;
b1a = 2.25*10^-12;
b2a = 22.475;
b0s = 3.9875*10^-4;
b1s = 2.225*10^-14;
b2s = 26.55;

TC = 1.1*max(T);
TC_it = 0;

while (TC-TC_it)^2>0.0000001
    
TC_it = TC;

% If using the model, uncomment these:
% Non-constant term

% erra = pea*exp(ea*T/TC);
% errs = pes*exp(es*T/TC);

% Constant term

erra = b0a + b1a*exp(b2a*T/TC);
errs = b0s + b1s*exp(b2s*T/TC);

errs = rhos * (1/beta) .* errs ./ delrho; % This should still be statistically accurate

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


end

% alpha = 0.95; % For normal use
alpha = 0.995; % For showing the nonphysical region

% CI = 1.3*tinv(alpha,n-4)/sqrt(n); % The 1.3 is to approximately make it the two-tailed solution in the range of 4-6 temperatures at 95%
CI = 1.42*tinv(alpha,n-4)/sqrt(n); % The 1.3 is to approximately make it the two-tailed solution for 3 temperatures at 99.5%

low_TC = TC-dTC*CI;
high_TC = TC+dTC*CI;
low_rhoc = rhoc-drhoc*CI;
high_rhoc = rhoc+drhoc*CI;

figure
hold
plot([low_rhoc,low_rhoc],[low_TC,high_TC])
plot([high_rhoc,high_rhoc],[low_TC,high_TC])
plot([low_rhoc,high_rhoc],[low_TC,low_TC])
plot([low_rhoc,high_rhoc],[high_TC,high_TC])
scatter(rhoc,TC,'b')
scatter(rhol,T,'g')
scatter(rhog,T,'g')
hold

