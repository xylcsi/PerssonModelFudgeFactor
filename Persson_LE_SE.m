function U_el = Persson_LE_SE(p_bar, E_star, ql, qr, xi, C0, H, option, misc)
% This function is a gateway function of multiple variants of Persson's theory
% for solving elastic strain energy. 
%
% Inputs: 
% p_bar  = Average contact pressure, [scalar/array];
% E_star = Plane strain modulus;
%   ql   = Lower cutoff wavenumber;
%   qr   = Roll-off wavenumber;
%   xi   = Scale (upper cutoff wavenumber = xi*ql);
%   C0   = Constant proportionality in the power-law form of PSD
%   H    = Hurst exponent;
% option = Model name;
%        = 'Persson01': No correction.
%        = 'YP08': Fudge factor provided by Yang and Persson (2008).
%        = 'WM17': Fudge factor provided by Wang and Muser (2017).
%        = 'Xu24': New fudge factor.
% misc   = miscellaneous parameters  [optional]
%        = [n_xi] if option = 'Persson01',
%        = [gamma; n_xi] if option = 'YP08', 'WM17, 'Xu24'.
%
% Outputs:
% U_el = Linear elastic strain energy [scalar]
%
% Syntax: 
% 
% U_el = Persson_LE_Area(p_bar, E_star, ql, qr, xi, C0, H, option)
% U_el = Persson_LE_Area(p_bar, E_star, ql, qr, xi, C0, H, option, gamma)
% U_el = Persson_LE_Area(p_bar, E_star, ql, qr, xi, C0, H, option, [gamma, n_xi])
%
% Debug log:
% 
% Notes: 
% 1. When option = 'WM17', gamma = 5/9.
% 2. When option = 'Xu24' or 'YP08', gamma = 0.45.
%
switch option
    case 'Persson01'
        if nargin == 8
            U_el = Persson_LE_SE_Persson01(p_bar, E_star, ql, qr, xi, C0, H);            
        elseif nargin == 9
            n_xi = misc(1); 
            U_el = Persson_LE_SE_Persson01(p_bar, E_star, ql, qr, xi, C0, H, n_xi);            
        end
    case 'YP08'
        if nargin == 8
            U_el = Persson_LE_SE_YP08(p_bar, E_star, ql, qr, xi, C0, H);
        else
            if max(size(misc)) == 1
                gamma = misc(1); 
                U_el = Persson_LE_SE_YP08(p_bar, E_star, ql, qr, xi, C0, H, gamma);
            elseif max(size(misc)) == 2
                gamma = misc(1); 
                n_xi = floor(misc(2)); 
                U_el = Persson_LE_SE_YP08(p_bar, E_star, ql, qr, xi, C0, H, gamma, n_xi);
            else
            end
        end
    case 'WM17'
        if nargin == 8
            U_el = Persson_LE_SE_WM17(p_bar, E_star, ql, qr, xi, C0, H);
        else
            if max(size(misc)) == 1
                gamma = misc(1); 
                U_el = Persson_LE_SE_WM17(p_bar, E_star, ql, qr, xi, C0, H, gamma);
            elseif max(size(misc)) == 2
                gamma = misc(1); 
                n_xi = floor(misc(2)); 
                U_el = Persson_LE_SE_WM17(p_bar, E_star, ql, qr, xi, C0, H, gamma, n_xi);
            else
            end
        end        
    case 'Xu24'
        if nargin == 8
            U_el = Persson_LE_SE_Xu24(p_bar, E_star, ql, qr, xi, C0, H);
        else
            if max(size(misc)) == 1
                gamma = misc(1); 
                U_el = Persson_LE_SE_Xu24(p_bar, E_star, ql, qr, xi, C0, H, gamma);
            elseif max(size(misc)) == 2
                gamma = misc(1); 
                n_xi = floor(misc(2)); 
                U_el = Persson_LE_SE_Xu24(p_bar, E_star, ql, qr, xi, C0, H, gamma, n_xi);
            else
            end
        end
    case 'Mix'
        if nargin == 8
            U_el = Persson_LE_SE_Mixture(p_bar, E_star, ql, qr, xi, C0, H);
        else
            if max(size(misc)) == 1
                gamma = misc(1); 
                U_el = Persson_LE_SE_Mixture(p_bar, E_star, ql, qr, xi, C0, H, gamma);
            elseif max(size(misc)) == 2
                gamma = misc(1); 
                n_xi = floor(misc(2)); 
                U_el = Persson_LE_SE_Mixture(p_bar, E_star, ql, qr, xi, C0, H, gamma, n_xi);
            else
            end
        end
    otherwise
        error('Wrong input of option! \n'); 
end