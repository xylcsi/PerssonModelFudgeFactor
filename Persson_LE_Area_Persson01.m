function [CR, Vp] = Persson_LE_Area_Persson01(p_bar, E_star, ql, qr, xi, C, H)
% This function calculates the relative contact area using the original
% Persson's theory [1, 2]. This function is based on the work of Persson (2001) 
% [1] (see also [2]) with no correction term.
% 
% Inputs: 
% p_bar  = Average contact pressure [scalar/array];
% E_star = Plane strain modulus;
%   ql   = Lower cutoff wavenumber;
%   qr   = Roll-off wavenumber;
%   xi   = Scale (upper cutoff wavenumber = xi*ql);
%   C    = Constant in the PSD;
%   H    = Hurst exponent.
%
% Outputs:
% CR = Relative contact area [scalar/array]
% Vp = Variance of contact pressure [scalar]
% 
% References:
% [1] Persson, B.N., 2001. Theory of rubber friction and contact mechanics. 
% The Journalof Chemical Physics,115(8), pp.3840-3861.
%
% [2] Manners, W. and Greenwood, J.A., 2006. Some observations on Persson’s 
% diffusiontheory of elastic contact. Wear, 261(5-6), pp.600-610.
% 
Vp = Variance_pressure(E_star, ql, qr, xi, C, H, 'V'); 
CR = erf(p_bar./sqrt(2*Vp)); % Original formula
