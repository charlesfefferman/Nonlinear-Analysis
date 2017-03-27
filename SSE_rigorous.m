function [ TC,rho_c ] = SSE_rigorous( TvL,pv,pL )
%This is the rigorous way to calculate SSE to be used in fminsearch

n = 2*length(TvL);

beta = 0.32;

rhoa = (pv + pL) / 2;

delrho = (pL - pv) / 2;

b0a = 5*10^-4;
b1a = 2.25*10^-12;
b2a = 22.475;
b0s = 3.9875*10^-4;
b1s = 2.225*10^-14;
b2s = 26.55;

TC_g = 1.1*max(TvL);
rhoc_g = 0.9*min(rhoa);
A_g = 4*10^-4;
B_g = 0.1;

P_g = [TC_g, rhoc_g, A_g, B_g];

erra = @(T,TC) b0a + b1a.*exp(b2a.*T/TC); % Error model for addition
errs = @(T,TC) b0s + b1s.*exp(b2s.*T/TC); % Error model for subtraction
LRD = @(T,TC,rho_c,A) rho_c + A.*(TC - T); % Law of rectilinear diameters
DSL = @(T,TC,B) B.*(TC-T).^beta; % Density scaling law

SSE = @(TC,rho_c,A,B) sum(((rhoa-LRD(TvL,TC,rho_c,A))./(erra(TvL,TC))).^2 + ((delrho - DSL(TvL,TC,B))./(errs(TvL,TC))).^2);

SSE_P = @(P) SSE(P(1),P(2),P(3),P(4));

% P = fminsearch(SSE_P,P_g);

P = fmincon(SSE_P,P_g,[],[],[],[],[1.01*max(TvL) 0.1*min(rhoa) 2*10^-4, 0],[]);

TC = P(1);
rho_c = P(2);


end

