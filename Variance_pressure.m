function res = Variance_pressure(E_star, ql, qr, xi, C0, H, item)
% Determine the varaince of the contact pressure distribution when the rough 
% surface is completely flattened.
%
% Inputs: 
% p_bar  = Average contact pressure [scalar/array];
% E_star = Plane strain modulus;
%   ql   = Lower cutoff wavenumber;
%   qr   = Roll-off wavenumber;
%   xi   = Scale (upper cutoff wavenumber = xi*ql);
%   C0   = Constant in the PSD;
%   H    = Hurst exponent.
%
% Outputs
% res   = Variance of contact pressure (item = 'V'); 
%       = Derivative of the variance of contact pressure w.r.t scale (item = 'dV/dxi'); 
%
if strcmp(item, 'V') % Variance
    if xi*ql < qr
        res = 0.125*pi*E_star^2*C0*qr^(-2*H-2)*ql^4*(xi^4 - 1); 
    else
        res = 0.125*pi*E_star^2*C0*(qr^(-2*H-2)*(qr^4 - ql^4) + ...
             2/(1 - H)*(xi^(-2*H+2)*ql^(-2*H+2) - qr^(-2*H+2))); 
    end
elseif strcmp(item, 'dV/dxi') % Derivative of varaince w.r.t scale
    res = pi/2*E_star^2*(xi*ql).^3*PSD_SASurf(ql, qr, xi*ql, xi*ql, C0, H)*ql; 
else
    error('Wrong inputs for Variance_pressure()!');
end