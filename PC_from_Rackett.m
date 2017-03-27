function [ PC ] = PC_from_Rackett( TC, rho_c, VC, A, B, TR)
%This function returns PC from the fits to the rectilinear and density
%scaling laws for a specific reduced temperature

beta = 0.32;

D = 2/7;
Rg = 8.314472;

rho_L = rho_c + A.*TC.*(1-TR) + B.*TC.^beta.*(1-TR).^beta;

PC = exp(log(rho_c./rho_L)./((1-TR).^D) - log(VC./Rg./TC));

% hold
% plot(TR*TC,rho_L,'r')
% hold

end

