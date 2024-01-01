function PSD = PSD_SASurf(ql, qr, qs, q, C0, H)
% This function calculates the PSD of a self-affine, isotropic, random, 
% bandwidth-limited rough surface. 
%
% Inputs: 
% ql = Lower cutoff wavenumber;
% qr = Roll-off wavenumber;
% qs = Upper cutoff wavenumber; 
% q  = Wavenumber at which the PSD is evaluated;
% C0 = Constant in the PSD;
% H  = Hurst exponent; 
% 
% Outputs: 
% PSD = Power spectrum density at the wavenumber q.
% 
if q >= ql && q <= qr
    PSD = C0*qr^(-2*(1 + H)); 
elseif q > qr && q <= qs
    PSD = C0*q^(-2*(1 + H)); 
else
    PSD = 0; 
end