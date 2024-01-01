function [CR, Vp] = Persson_LE_Area_WM17(p_bar, E_star, ql, qr, xi, C0, H, gamma, n_xi)
% This function calculates the relative contact area using Persson's theory with
% the fudge factor proposed by Wang and Muser [1]. 
% 
% Inputs: 
% p_bar  = Average contact pressure [scalar/array];
% E_star = Plane strain modulus;
%   ql   = Lower cutoff wavenumber;
%   qr   = Roll-off wavenumber;
%   xi   = Scale (upper cutoff wavenumber = xi*ql);
%   C    = Constant in the PSD;
%   H    = Hurst exponent;
%  gamma = Fitting parameter; [optional]
%  n_xi  = Number of grid on the scale axis. [optional]
%
% Outputs:
% CR = Relative contact area. [scalar/array]
% 
% Syntax:
% 
%    [CR, Vp] = Persson_LE_Area_WM17(p_bar, E_star, fl, fr, xi, C, H)
%    [CR, Vp] = Persson_LE_Area_WM17(p_bar, E_star, fl, fr, xi, C, H, gamma)
%    [CR, Vp] = Persson_LE_Area_WM17(p_bar, E_star, fl, fr, xi, C, H, gamma, n_xi)
%
% References:
% [1] Wang A, Muser M H. Gauging Persson theory on adhesion[J]. 
%     Tribology Letters, 2017, 65(3): 1 – 10.
% 
if nargin == 7
    n_xi = 500;
    gamma = 5/9; 
end
if nargin == 8
    n_xi = 500; 
end
       n = max(size(p_bar));
      CR = zeros(n, 1); 
      Vp = zeros(n, 1); % Vp at partial contact
     Vpc = zeros(n_xi, 1); % Vp at complete contact
 C_array = zeros(n_xi, 1); % PSD of rough surface
xi_array = logspace(0, log10(xi), n_xi); 
xi_array = xi_array(:); 
q_array  = xi_array*ql; 
Ar_array = zeros(n_xi, 1); 
Ar_array(1) = 1;
for i = 1: n_xi
    C_array(i) = PSD_SASurf(ql, qr, xi*ql, q_array(i), C0, H); 
        Vpc(i) = Variance_pressure(E_star, ql, qr, xi_array(i), C0, H, 'V'); 
end
%
func = @(Ar, a, b, c, p_bar) Ar - erf(p_bar./sqrt(a + b*Ar.^2 + c*Ar.^4));  
for m = 1: n
    if p_bar(m) == 0
        CR(m) = 0; 
    else
        for i = 2: n_xi            
            func1 = @(Ar, a, b, c) func(Ar, a, b, c, p_bar(m)); 
            a = 2*Vp(m) + 2*gamma*(Vpc(i) - Vpc(i-1)); 
            b = -4/9*(Vpc(i) - Vpc(i-1)); 
            c = 2*(11/9 - gamma)*(Vpc(i) - Vpc(i-1)); 
            func2 = @(Ar) func1(Ar, a, b, c); 
            Ar_array(i) = fzero(func2, Ar_array(i-1)); 
            Vp(m) = 0.5*(a + b*Ar_array(i)^2 + c*Ar_array(i)^4); 
        end
        CR(m) = Ar_array(end); 
        Ar_array = zeros(n_xi, 1); 
        Ar_array(1) = 1; 
    end
end