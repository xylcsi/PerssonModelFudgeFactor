function Uel = Persson_LE_SE_YP08(p_bar, E_star, ql, qr, xi, C0, H, gamma, n_xi)
% This function calculates the strain energy using Persson's theory with the
% fudge factor provided by Yang and Persson [1-3].  
% 
% Inputs: 
% p_bar  = Average contact pressure [scalar];
% E_star = Plane strain modulus;
%   ql   = Lower cutoff wavenumber;
%   qr   = Roll-off wavenumber;
%   xi   = Scale (upper cutoff wavenumber = xi*ql);
%   C0   = Constant in the PSD;
%   H    = Hurst exponent;
%  gamma = Fitting parameter; [optional]
%  n_xi  = Number of grid on the scale axis. [optional]
%
% Outputs:
% SE = Linear elastic strain per unit area. [scalar]
% 
% References:
% [1] Yang C, Persson B. Contact mechanics: contact area and interfacial 
% separationfrom small contact to full contact[J]. Journal of Physics: Condensed 
% Matter, 2008,20(21) : 215214.
%
% [2] Almqvist, A., Campana, C., Prodanov, N. and Persson, B.N.J., 2011. 
% Interfacialseparation between elastic solids with randomly rough surfaces: 
% comparison between theory and numerical techniques. Journal of the Mechanics 
% and Physics ofSolids, 59(11), pp.2355-2369.
% 
% [3] Afferrante L, Bottiglione F, Putignano C, et al. Elastic contact mechanics 
% of randomly rough surfaces: an assessment of advanced asperity models and 
% Persson%stheory[J]. Tribology Letters, 2018, 66(2) : 1 – 13.
%
if nargin == 7
    n_xi = 500;
    gamma = 0.42; 
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
func = @(Ar, a, b) Ar - erf(p_bar./sqrt(a + b*Ar.^2));  

if p_bar == 0
    Uel = 0; 
else
    Vp  = 0; 
    Uel = (q_array(1))^2*C_array(1)*(q_array(2) - q_array(1)); 
    for i = 2: n_xi            
        a = 2*Vp + 2*gamma*(Vpc(i) - Vpc(i-1)); 
        b = 2*(1 - gamma)*(Vpc(i) - Vpc(i-1)); 
        func1 = @(Ar) func(Ar, a, b); 
        Ar_array(i) = fzero(func1, Ar_array(i-1)); 
        S = gamma + (1 - gamma)*Ar_array(i)^2; 
        Uel = Uel + (q_array(i))^2*C_array(i)*Ar_array(i)*S*(q_array(i) - q_array(i - 1)); 
        Vp = 0.5*(a + b*Ar_array(i)^2); 
    end
    Uel = pi*E_star/2*Uel; 
end
