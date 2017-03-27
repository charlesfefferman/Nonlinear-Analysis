function [ rho_L ] = rho_L_hat( T, rhoc, TC, A, B, beta )
% Uses rectilinear and scaling to predict LDN

rho_L = rhoc + A.*(TC - T) + B.*(TC - T).^beta;


end

