function Uel = Persson_LE_SE_Persson01(p_bar, E_star, ql, qr, xi, C0, H, n_xi)
% This function calculates the strain energy using uncorrected Persson's theory
% given by Persson [1, 2]
% 
% Inputs: 
% p_bar  = Average contact pressure [scalar];
% E_star = Plane strain modulus;
%   ql   = Lower cutoff wavenumber;
%   qr   = Roll-off wavenumber;
%   xi   = Scale (upper cutoff wavenumber = xi*ql);
%   C0   = Constant in the power-law form of PSD;
%   H    = Hurst exponent;
%  n_xi  = Number of grid on the scale axis. [optional]
%
% Outputs:
% SE = Linear strain energy per unit area [scalar].
% 
% References:
% [1] Persson B. Relation between interfacial separation and load: a general 
% theory of contact mechanics[J]. Physical review letters, 2007, 99(12) : 125502.
% 
% [2] Persson B. On the elastic energy and stress correlation in the contact 
% between elastic solids with randomly rough surfaces[J]. Journal of Physics: 
% Condensed Matter, 2008, 20(31) : 312001.
%
if nargin == 7
    n_xi = 500; 
end
xi_array = logspace(0, log10(xi), n_xi); 
PSD = PSD_SASurf(ql, qr, xi*ql, xi_array(1)*ql, C0, H); 
CR = Persson_LE_Area_Persson01(p_bar, E_star, ql, qr, xi_array(1)*ql, C0, H); 
Uel = (xi_array(1)*ql)^2*PSD*CR*(xi_array(2) - xi_array(1))*ql; 
for i = 2: n_xi
     CR = Persson_LE_Area_Persson01(p_bar, E_star, ql, qr, xi_array(i)*ql, C0, H); 
    PSD = PSD_SASurf(ql, qr, xi*ql, xi_array(i)*ql, C0, H);
    Uel = Uel + (xi_array(i)*ql)^2*PSD*CR*(xi_array(i) - xi_array(i-1))*ql; 
end
Uel = pi*E_star/2*Uel; 