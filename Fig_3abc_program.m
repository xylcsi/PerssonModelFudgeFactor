clc
clear all
close all
% 
% This script recreates Fig. 3(a)-(c) in the present work. 
% 
% Problem statement
E_star = 1; % (Pa) Plane strain modulus
ql = 1e4; % (1/m) Lower cut-off frequency
qr = ql; % (1/m) Roll-off frequency
qs = ql*100; % (1/m) Upper cut-off frequency
xi = qs/ql; % Scale
H  = 0.8; % Hurst dimension
h_rms = 6e-6; % (m) root mean square roughness
C0 = h_rms^2*H/pi/(ql^(-2*H) - qs^(-2*H)); % Constant proportionality of PSD
%
% 2. PDF of contact pressure
p_bar = 0.1130; % Pa
%
[p_Persson01, PDF_pres_Persson01] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'Persson01', [200; 6]);
[p_YP08, PDF_pres_YP08] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'YP08', [200; 6; 0.45; 500]);
[p_WM17, PDF_pres_WM17] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'WM17', [200; 6; 5/9; 500]);
[p_Xu24, PDF_pres_Xu24] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'Xu24', [200; 6; 0.45]);
load('Data/PDF0_113_GFMD_Zhou.mat'); 
%
figure; 
hold on
plot(p_Persson01, PDF_pres_Persson01, '-', 'LineWidth', 2);
plot(p_YP08, PDF_pres_YP08, '-', 'LineWidth', 2);
plot(p_WM17, PDF_pres_WM17, '-', 'LineWidth', 2);
plot(p_Xu24, PDF_pres_Xu24, '-', 'LineWidth', 2); 
errorbar(p_bin, PDF0_113_mean, PDF0_113_std, 'o');
hold off
xlabel('$p$ (Pa)', 'interpreter', 'latex'); 
ylabel('$P_0(p, \xi)$', 'interpreter', 'latex'); 
legend('Persson, 2001', 'Yang and Persson, 2008', 'Wang and Muser, 2017', 'Present work', 'GFMD'); 

% save('Fig_R2_pbar0_113.mat', 'p_Persson01', 'PDF_pres_Persson01', 'p_YP08', 'PDF_pres_YP08', ...
%      'p_WM17', 'PDF_pres_WM17', 'p_Xu24', 'PDF_pres_Xu24', 'p_bin', 'PDF0_113_mean', 'PDF0_113_std'); 
%
p_bar = 0.24; % Pa
%1
[p_Persson01, PDF_pres_Persson01] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'Persson01', [200; 6]);
[p_YP08, PDF_pres_YP08] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'YP08', [200; 6; 0.45; 500]);
[p_WM17, PDF_pres_WM17] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'WM17', [200; 6; 5/9; 500]);
[p_Xu24, PDF_pres_Xu24] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'Xu24', [200; 6; 0.45]);
load('Data/PDF0_24_GFMD_Zhou.mat'); 
%
figure; 
hold on
plot(p_Persson01, PDF_pres_Persson01, '-', 'LineWidth', 2);
plot(p_YP08, PDF_pres_YP08, '-', 'LineWidth', 2);
plot(p_WM17, PDF_pres_WM17, '-', 'LineWidth', 2);
plot(p_Xu24, PDF_pres_Xu24, '-', 'LineWidth', 2); 
errorbar(p_bin, PDF0_24_mean, PDF0_24_std, 'o');
hold off
xlabel('$p$ (MPa)', 'interpreter', 'latex'); 
ylabel('$P_0(p, \xi)$', 'interpreter', 'latex'); 
legend('Persson, 2001', 'Yang and Persson, 2008', 'Wang and Muser, 2017', 'Present work', 'GFMD'); 
% save('Fig_R2_pbar0_24.mat', 'p_Persson01', 'PDF_pres_Persson01', 'p_YP08', 'PDF_pres_YP08', ...
%      'p_WM17', 'PDF_pres_WM17', 'p_Xu24', 'PDF_pres_Xu24', 'p_bin', 'PDF0_24_mean', 'PDF0_24_std'); 
%
p_bar = 0.025; % Pa
%
[p_Persson01, PDF_pres_Persson01] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'Persson01', [200; 6]);
[p_YP08, PDF_pres_YP08] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'YP08', [200; 6; 0.45; 500]);
[p_WM17, PDF_pres_WM17] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'WM17', [200; 6; 5/9; 500]);
[p_Xu24, PDF_pres_Xu24] = Persson_LE_PDF_Pres(p_bar, E_star, ql, qr, xi, C0, H, 'Xu24', [200; 6; 0.45]);
load('Data/PDF0_025_GFMD_Zhou.mat'); 
%
figure; 
hold on
plot(p_Persson01, PDF_pres_Persson01, '-', 'LineWidth', 2);
plot(p_YP08, PDF_pres_YP08, '-', 'LineWidth', 2);
plot(p_WM17, PDF_pres_WM17, '-', 'LineWidth', 2);
plot(p_Xu24, PDF_pres_Xu24, '-', 'LineWidth', 2); 
errorbar(p_bin, PDF0_025_mean, PDF0_025_std, 'o');
hold off
xlabel('$p$ (Pa)', 'interpreter', 'latex'); 
ylabel('$P_0(p, \xi)$', 'interpreter', 'latex'); 
legend('Persson, 2001', 'Yang and Persson, 2008', 'Wang and Muser, 2017', 'Present work', 'GFMD');  
% save('Fig_R2_pbar0_025.mat', 'p_Persson01', 'PDF_pres_Persson01', 'p_YP08', 'PDF_pres_YP08', ...
%      'p_WM17', 'PDF_pres_WM17', 'p_Xu24', 'PDF_pres_Xu24', 'p_bin', 'PDF0_025_mean', 'PDF0_025_std'); 