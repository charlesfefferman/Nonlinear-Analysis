clear

% In this code, the "traditional" in the title is referring to the
% traditional amount of time letting the simulation run and the rigorous
% propagation means that the rectilinear and scaling densities were
% followed independently.  This script can not be used for literature
% results because this requires the raw data.

% This is my TraPPE octane data
% Towhee gave:
% rhoc 0.239147 0.02396
% TC 569.312433 24.019290

TC = 569;

T = [390 440 490 515 543];
rhog = [0.004219 0.013343 0.032864 0.05073475 0.0840558];
rhol = [0.625876 0.573911 0.511936 0.475198579 0.418986443];
errg = [0.000625 0.002313 0.005560 0.008574752 0.013940418];
errl = [0.003440 0.003881 0.004561 0.009607261 0.018147998];
erra = [0.001007704 0.0024603 0.0043773 0.0086883 0.013896621];	% Found using traditional length of simulation
errs = [0.001047761 0.0020383 0.0025878 0.0027253 0.008290142];

beta = 0.32;

rhoa = (rhog + rhol) / 2;

delrho = (rhol - rhog) / 2;
rhos = delrho;
rhos = rhos.^(1/beta);

TC_it = 0;

while (TC-TC_it)^2>0.0000001
    
TC_it = TC;

errs_alt = rhos * (1/beta) .* errs ./ delrho; % This should still be statistically accurate

[as, bs, sigas, sigbs] = towhee_fit(T,rhos,errs_alt);

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