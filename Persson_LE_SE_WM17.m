function Uel = Persson_LE_SE_WM17(p_bar, E_star, ql, qr, xi, C0, H, gamma, n_xi)
% This function calculates the strain energy using Persson's theory with the
% fudge factor provided by Wang and Muser [1].  
% 
% Inputs: 
% p_bar  = Average contact pressure [scalar];
% E_star = Plane strain modulus;
%   ql   = Lower cutoff wavenumber;
%   qr   = Roll-off wavenumber;
%   xi   = Scale (upper cutoff wavenumber = xi*ql);
%   C0   = Constant in the power-law form of PSD;
%   H    = Hurst exponent.
%  gamma = Fitting parameter. [optional]
%  n_xi  = Number of grid on the scale axis [optional]
%
% Outputs:
% SE = Linear elastic strain per unit area [scalar]
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
%
     Vpc = zeros(n_xi, 1); % Pressure variance at complete contact
 C_array = zeros(n_xi, 1); % PSD of rough surface
Ar_array = zeros(n_xi, 1); % Relative contact area array
Ar_array(1) = 1; 
xi_array = logspace(0, log10(xi), n_xi); 
xi_array = xi_array(:); 
q_array  = xi_array*ql; 
for i = 1: n_xi
    C_array(i) = PSD_SASurf(ql, qr, xi*ql, q_array(i), C0, H); 
        Vpc(i) = Variance_pressure(E_star, ql, qr, xi_array(i), C0, H, 'V'); 
end
%
func = @(Ar, a, b, c) Ar - erf(p_bar./sqrt(a + b*Ar.^2 + c*Ar.^4));  

if p_bar == 0
    Uel = 0; 
else
    Uel = (q_array(1))^2*C_array(1)*(q_array(2) - q_array(1)); 
     Vp = 0; 
    for i = 2: n_xi            
        a = 2*Vp + 2*gamma*(Vpc(i) - Vpc(i-1)); 
        b = -4/9*(Vpc(i) - Vpc(i-1)); 
        c = 2*(11/9 - gamma)*(Vpc(i) - Vpc(i-1)); 
        func1 = @(Ar) func(Ar, a, b, c); 
        Ar_array(i) = fzero(func1, Ar_array(i-1)); 
        S = gamma - 2/9*Ar_array(i)^2 + (11/9 - gamma)*Ar_array(i)^4; 
        Uel = Uel + (q_array(i))^2*C_array(i)*Ar_array(i)*S*(q_array(i) - q_array(i - 1)); 
        Vp = 0.5*(a + b*Ar_array(i)^2 + c*Ar_array(i)^4);
    end
    Uel = pi*E_star/2*Uel; 
end
