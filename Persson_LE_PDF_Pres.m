function [p, PDF, Vp] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, option, misc)
% This function is a gateway function of multiple variants of Persson's theory
% for solving PDF of contact pressure. 
%
% Inputs: 
% p_bar  = Average contact pressure, [scalar]
% E_star = Plane strain modulus,
%   ql   = Lower cut-off wavenumber;
%   qr   = Roll-off wavenumber;
%  xi    = Scale (upper cut-off wavenumber = xi*ql);
%  C0    = Constant in the PSD;
%  H     = Hurst exponent;
% option = Model name;
%        = 'Persson01': No correction.
%        = 'YP08': Fudge factor provided by Yang and Persson (2008).
%        = 'WM17': Fudge factor provided by Wang and Muser (2017).
%        = 'Xu24': New fudge factor.
% misc   = miscellaneous parameters  [optional]
%        = [np; a] if option = 'Persson01'
%        = [np; a; gamma; n_xi] if option = 'YP08'
%        = [np; a; gamma; n_xi] if option = 'WM17'
%        = [np; a; gamma] if option = 'Xu23'
%
% Outputs:
% p   = Contact pressure [array]
% PDF = PDF of contact pressure [array]
%
% Syntax: 
% 
%    [p, PDF] = Persson_LE_PDF_pres(p_bar, E_star, ql, qr, xi, C0, H, option)
%    [p, PDF] = Persson_LE_PDF_pres(p_bar, E_star, ql, qr, xi, C0, H, option, misc)
%
% Debug log:
% 
switch option
    case 'Persson01'
        Vp = Variance_pressure(E_star, ql, qr, xi, C0, H, 'V'); 
        if nargin == 8
            np = 200; 
             a = 6; 
        else
            np = misc(1); 
            a  = misc(2); 
        end
    case 'YP08'
        if nargin == 8
            [~, Vp] = Persson_LE_Area(p_bar, E_star, ql, qr, xi, C0, H, 'YP08'); 
            np = 200; 
             a = 6; 
        else % nargin == 9    
            [~, Vp] = Persson_LE_Area(p_bar, E_star, ql, qr, xi, C0, H, 'YP08', misc(3:4)); 
            np = misc(1); 
            a  = misc(2); 
        end
    case 'WM17'
        if nargin == 8
            [~, Vp] = Persson_LE_Area(p_bar, E_star, ql, qr, xi, C0, H, 'WM17'); 
            np = 200; 
             a = 6; 
        else % nargin == 9
            [~, Vp] = Persson_LE_Area(p_bar, E_star, ql, qr, xi, C0, H, 'WM17', misc(3:4)); 
            np = misc(1); 
            a  = misc(2); 
        end
    case 'Xu24'
        if nargin == 8
            [~, Vp] = Persson_LE_Area(p_bar, E_star, ql, qr, xi, C0, H, 'Xu24'); 
            np = 200; 
             a = 6; 
        else
            [~, Vp] = Persson_LE_Area(p_bar, E_star, ql, qr, xi, C0, H, 'Xu24', misc(3)); 
            np = misc(1); 
            a  = misc(2); 
        end
    otherwise 
        error('Wrong input of option! \n'); 
end
  p = linspace(0, p_bar + a*sqrt(Vp), np); 
PDF = 1/sqrt(2*pi*Vp)*(exp(-(p - p_bar).^2/2/Vp) - exp(-(p + p_bar).^2/2/Vp));   
end