clc
clear
close all

% This script recreates Fig. 3d in the present work. 
data = load('Data/Wang_Muser_Fig2_data.txt'); 
load_GFMD_WM17 = [0; data(:, 1)]; 
PDF_pres_GFMD_WM17 = [0; data(:, 2)]; 

% Surface PSD
    E  = 25; % (MPa) 
    nu = 0; % 
E_star = E/(1 - nu^2);
    Lx = 0.1; % (mm)
    Ly = 0.1; % (mm)
    ql = 2*pi/Lx; % Lower cutoff wavenumber (1/mm)
    qr = 2*pi/0.020; % Roll-off wavenumber (1/mm)
    qs = 1/1e-4; % Upper cutoff wavenumber (1/mm)
    xi = qs/ql; % Scale     
     H = 0.8; % Hurst dimension
 g_bar = 1; % R.M.S gradient of the surface
    C0 = g_bar*2/pi/(qr^(-2*(1 + H))* ...
         (qr^4 - ql^4) + 2/(1 - H)*(qs^(2 - 2*H) - qr^(2 - 2*H))); 
 p_bar = 0.25; % MPa
%
[p_Persson01, PDF_pres_Persson01] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'Persson01');
[p_YP08, PDF_pres_YP08] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'YP08', [200; 6; 0.45; 500]);
[p_WM17, PDF_pres_WM17] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'WM17', [200; 6; 5/9; 500]);
[p_Xu24, PDF_pres_Xu24] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'Xu24', [200; 6; 0.45]);
%
figure; 
hold on
plot(p_Persson01/E_star, PDF_pres_Persson01*E_star, 'ro', 'LineWidth', 2);
plot(p_YP08/E_star, PDF_pres_YP08*E_star, 'gd', 'LineWidth', 2);
plot(p_WM17/E_star, PDF_pres_WM17*E_star, 'bx', 'LineWidth', 2);
plot(p_Xu24/E_star, PDF_pres_Xu24*E_star, 'k-', 'LineWidth', 2); 
plot(load_GFMD_WM17, PDF_pres_GFMD_WM17, 'ko', 'LineWidth', 2);
hold off
xlabel('$p/E^*$', 'interpreter', 'latex'); 
ylabel('$P_0(p/E^*, \xi)$', 'interpreter', 'latex'); 
xlim([0, 1.3]);
ylim([0, 0.045]); 
legend('Persson, 2001', ...
       'Yang and Persson, 2008', ...
       'Wang and Muser, 2017', 'Present work', 'GFMD'); 
p_Persson01 = p_Persson01/E_star; 
PDF_pres_Persson01 = PDF_pres_Persson01*E_star; 
p_YP08 = p_YP08/E_star; 
PDF_pres_YP08 = PDF_pres_YP08*E_star; 
p_WM17 = p_WM17/E_star; 
PDF_pres_WM17 = PDF_pres_WM17*E_star; 
p_Xu24 = p_Xu24/E_star; 
PDF_pres_Xu24 = PDF_pres_Xu24*E_star; 
p_GFMD_WM17 = load_GFMD_WM17; 
% save('Fig_R2_d.mat', 'p_Persson01', 'PDF_pres_Persson01', 'p_YP08', 'PDF_pres_YP08', ...
%      'p_WM17', 'PDF_pres_WM17', 'p_Xu24', 'PDF_pres_Xu24', 'p_GFMD_WM17', 'PDF_pres_GFMD_WM17');

