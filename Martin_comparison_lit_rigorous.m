clear

% Ethane

T = [178 217 236 256];

% AMBER

pv = [0.0032 0.0179 0.0291 0.054];
pL = [0.528 0.466 0.424 0.364];

% CHARMM

% pv = [0.0017 0.0104 0.0214 0.0441];
% pL = [0.558 0.502 0.468 0.430];
% 

% OPLS

pv = [0.0025 0.0143 0.0263 0.057];
pL = [0.551 0.487 0.447 0.389];

% % TraPPE
% 
% pv = [0.0021 0.0111 0.0204 0.0350];
% pL = [0.551 0.499 0.47 0.434];
% 
% % Validation data
% 
T = [178 197 217 256 275 279 283 288];
pv = [0.0023 0.0053 0.0111 0.0350 0.0598 0.0648 0.0739 0.090];
pL = [0.5512 0.5262 0.4984 0.4342 0.3937 0.3835 0.3726 0.3589];
T = T(5:end);
pv = pv(5:end);
pL = pL(5:end);
T = [T T T T T T T T];
pv = [pv pv pv pv pv pv pv pv];
pL = [pL pL pL pL pL pL pL pL];

[TC,pc] = towhee_error_model(T,pv,pL);
    
% [TC_low,TC_high,pc_low,pc_high] = rigorous_statistics(T,pv,pL)


