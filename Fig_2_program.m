clc
clear all
close all
% 
% This script recreates Fig. 2 in the present work. 
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
% 1. Relative contact ratio
p_bar_max = 5e-1;
p_bar_min = 5e-4; 
np = 300; 
p_bar_dense = logspace(log10(p_bar_min), log10(p_bar_max), np)'; 
%
CR_Persson01 = Persson_LE_Area(p_bar_dense, E_star, ql, qr, xi, C0, H, 'Persson01'); 
CR_YP08 = Persson_LE_Area(p_bar_dense, E_star, ql, qr, xi, C0, H, 'YP08', [0.45, 500]); 
CR_WM17 = Persson_LE_Area(p_bar_dense, E_star, ql, qr, xi, C0, H, 'WM17', [5/9, 500]); 
CR_Xu24 = Persson_LE_Area(p_bar_dense, E_star, ql, qr, xi, C0, H, 'Xu24', 0.45);
load('Data/Area_GFMD_Zhou.mat');
%
figure; 
hold on
plot([0; p_bar_dense], [0; CR_Persson01], '-', 'LineWidth', 2); 
plot([0; p_bar_dense], [0; CR_YP08], '-', 'LineWidth', 2); 
plot([0; p_bar_dense], [0; CR_WM17], '-', 'LineWidth', 2);
plot([0; p_bar_dense], [0; CR_Xu24], '-', 'LineWidth', 2); 
errorbar(p_bar, area_mean, area_std, 'o', 'LineWidth', 2);
hold off
xlabel('$\bar{p}$ (Pa)', 'Interpreter', 'latex'); 
ylabel('$A^*$', 'Interpreter', 'latex'); 
% set(gca,'Xscale','log');
% set(gca,'Yscale','log');
legend('Persson01', 'YP08', 'WM17', 'Present work', 'GFMD');
p_bar_dense = [0; p_bar_dense];
CR_Persson01 = [0; CR_Persson01]; 
CR_YP08 = [0; CR_YP08];
CR_WM17 = [0; CR_WM17]; 
CR_Xu24 = [0; CR_Xu24]; 
%save('Fig_R1.mat', 'p_bar_dense', 'CR_Persson01', 'CR_YP08', 'CR_WM17', 'CR_Xu24', 'p_bar', 'area_mean', 'area_std'); 