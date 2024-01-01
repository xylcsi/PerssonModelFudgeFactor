function [CR, Vp] = Persson_LE_Area(p_bar, E_star, ql, qr, xi, C0, H, option, misc)
% This function is a gateway function of multiple variants of Persson's theory
% for solving the relative contact area. 
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
%        = [] if option = 'Persson01',
%        = [gamma; n_xi] if option = 'YP08' or 'WM17,
%        = [gamma] if option = 'Xu24'.
%
% Outputs:
% CR = Contact ratio [scalar/array]
% Vp = Variance of contact pressure
%
% Syntax: 
% 
% CR = Persson_LE_Area(p_bar, E_star, ql, qr, xi, C0, H, option)
% CR = Persson_LE_Area(p_bar, E_star, ql, qr, xi, C0, H, option, [gamma])
% CR = Persson_LE_Area(p_bar, E_star, ql, qr, xi, C0, H, option, [gamma, n_xi])
%
% Debug log:
% 
% Notes: 
% 1. When option = 'WM17', gamma = 5/9.
% 2. When option = 'Xu24' or 'YP08', gamma = 0.45.
%
switch option
    case 'Persson01'
        [CR, Vp] = Persson_LE_Area_Persson01(p_bar, E_star, ql, qr, xi, C0, H);
    case 'YP08'
        if nargin == 8
            [CR, Vp] = Persson_LE_Area_YP08(p_bar, E_star, ql, qr, xi, C0, H);
        else
            if max(size(misc)) == 1
                gamma = misc(1); 
                [CR, Vp] = Persson_LE_Area_YP08(p_bar, E_star, ql, qr, xi, C0, H, gamma);
            elseif max(size(misc)) == 2
                gamma = misc(1); 
                n_xi = floor(misc(2)); 
                [CR, Vp] = Persson_LE_Area_YP08(p_bar, E_star, ql, qr, xi, C0, H, gamma, n_xi);
            else
            end
        end
    case 'WM17'
        if nargin == 8
            [CR, Vp] = Persson_LE_Area_WM17(p_bar, E_star, ql, qr, xi, C0, H);
        else
            if max(size(misc)) == 1
                gamma = misc(1); 
                [CR, Vp] = Persson_LE_Area_WM17(p_bar, E_star, ql, qr, xi, C0, H, gamma);
            elseif max(size(misc)) == 2
                gamma = misc(1); 
                n_xi = floor(misc(2)); 
                [CR, Vp] = Persson_LE_Area_WM17(p_bar, E_star, ql, qr, xi, C0, H, gamma, n_xi);
            else
            end
        end        
    case 'Xu24'
        if nargin == 8
            [CR, Vp] = Persson_LE_Area_Xu24(p_bar, E_star, ql, qr, xi, C0, H);
        else
            gamma = misc(1); 
            [CR, Vp] = Persson_LE_Area_Xu24(p_bar, E_star, ql, qr, xi, C0, H, gamma);
        end
    otherwise
        error('Wrong input of option! \n'); 
end