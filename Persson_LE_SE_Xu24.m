function Uel = Persson_LE_SE_Xu24(p_bar, E_star, ql, qr, xi, C0, H, gamma, n_xi)
% This function calculates the linear elastic strain energy stored in the 
% deformed interface under the framework of Persson's theory. 
% 
% Inputs: 
% p_bar  = Average contact pressure [scalar];
% E_star = Plane strain modulus;
%   ql   = Lower cutoff wavenumber;
%   qr   = Roll-off frequency;
%   xi   = Scale (upper cutoff wavenumber = xi*ql);
%   C    = Constant in the PSD;
%   H    = Hurst exponent.
%  gamma = Fitting parameter. [optional]
%  n_xi  = Number of grid on the scale axis [optional]
%
% Outputs:
% SE = Linear elastic strain per unit area [scalar]
% 
if nargin == 7
    n_xi = 500;
    gamma = 0.42; 
end
if nargin == 8
    n_xi = 500; 
end
%
if p_bar == 0
    Uel = 0; 
else
    Vpc_array = zeros(n_xi, 1); % Pressure variance at complete contact
     C_array = zeros(n_xi, 1); % PSD of rough surface
    CR_array = zeros(n_xi, 1); % Relative contact area array
    CR_array(1) = 1; 
    xi_array = logspace(0, log10(xi), n_xi); 
    xi_array = xi_array(:); 
    q_array  = xi_array*ql; 
    for i = 1: n_xi
        C_array(i) = PSD_SASurf(ql, qr, xi*ql, q_array(i), C0, H); 
      Vpc_array(i) = Variance_pressure(E_star, ql, qr, xi_array(i), C0, H, 'V'); 
       CR_array(i) = Persson_LE_Area_Xu24(p_bar, E_star, ql, qr, xi_array(i), C0, H, gamma); 
    end
    % Uel: part I
    Uel1 = q_array(1)^2*C_array(1)*(q_array(2) - q_array(1)); 
    for i = 2: n_xi
        S = gamma + (1 - gamma)*CR_array(i); 
        Uel1 = Uel1 + q_array(i)^2*C_array(i)*CR_array(i)*S*(q_array(i) - q_array(i-1)); 
    end
    Uel1 = pi*E_star/2*Uel1;
    % Uel: part II
    Uel2 = 0; 
    for i = 2: n_xi
        Uel2 = Uel2 + 1/q_array(i)*CR_array(i)*Vpc_array(i)*(CR_array(i) - CR_array(i-1)); 
    end
    Uel2 = (1 - gamma)/E_star*Uel2;
    Uel = Uel1 + Uel2; 
end