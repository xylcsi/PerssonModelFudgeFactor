function [CR, Vp] = Persson_LE_Area_Xu24(p_bar, E_star, ql, qr, xi, C, H, gamma)
% This function calculates the relative contact area using Persson's theory with
% the new fudge factor proposed in the present study.
% 
% Inputs: 
% p_bar  = Average contact pressure [scalar/array];
% E_star = Plane strain modulus;
%   ql   = Lower cutoff wavenumber;
%   qr   = Roll-off wavenumber;
%   xi   = Scale (upper cutoff wavenumber = xi*ql);
%   C    = Constant in the PSD;
%   H    = Hurst exponent;
% gamma  = Fitting parameter. [optional].
%
% Outputs:
% CR = Contact ratio [scalar/array].
% Vp = Variance of contact pressure [scalar/array]
% 
% Syntax:
% 
%   CR = Persson_LE_Area_Xu23(p_bar, E_star, fl, fr, xi, C, H)
%   CR = Persson_LE_Area_Xu23(p_bar, E_star, fl, fr, xi, C, H, gamma)
%
if nargin == 7
    gamma = 0.42; 
end
  Vpc = Variance_pressure(E_star, ql, qr, xi, C, H, 'V'); 
func = @(Ar, p_bar) Ar - erf(p_bar/sqrt(2*(gamma + (1 - gamma)*Ar)*Vpc));  
   n = max(size(p_bar)); 
  CR = zeros(n, 1); 
  Vp = zeros(n, 1); 
for i = 1: n
    func1 = @(Ar) func(Ar, p_bar(i)); 
    CR(i) = fzero(func1, erf(p_bar(i)/sqrt(2*Vpc)));
    Vp(i) = (gamma + (1 - gamma)*CR(i))*Vpc; 
end
