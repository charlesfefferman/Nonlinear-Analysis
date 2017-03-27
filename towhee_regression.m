function [ TC,rhoc,dTC,drhoc,A,b ] = towhee_regression(T,rhog,rhol,errg,errl )
%Towhee's algorithm

beta = 0.325;

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

A = b1/TC;
b = a0/(TC^beta);
end

