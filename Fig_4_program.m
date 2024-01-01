clc
clear all
close all
% 
% This script recreates Fig. 4 in the present work. 
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
% Linear elastic strain energy
%
p_bar_max = 5e-1;
p_bar_min = 5e-6; 
np = 300; 
p_bar_dense = logspace(log10(p_bar_min), log10(p_bar_max), np)'; 
for i = 1: np
    SE_Persson01(i) = Persson_LE_SE(p_bar_dense(i), E_star, ql, qr, xi, C0, H, 'Persson01'); 
         SE_YP08(i) = Persson_LE_SE(p_bar_dense(i), E_star, ql, qr, xi, C0, H, 'YP08'); 
         SE_WM17(i) = Persson_LE_SE(p_bar_dense(i), E_star, ql, qr, xi, C0, H, 'WM17');
         SE_Xu24(i) = Persson_LE_SE(p_bar_dense(i), E_star, ql, qr, xi, C0, H, 'Xu24');
          SE_Mix(i) = Persson_LE_SE(p_bar_dense(i), E_star, ql, qr, xi, C0, H, 'Mix');
end
load('Data/Uel_GFMD_Zhou.mat'); 
%
% Sepcial solution: Uel at complete contact
SE_CC = C0*qr^(-2*(1 + H))/3*(qr^3 - ql^3) + C0/(-2*H + 1)*(qs^(-2*H + 1) - qr^(-2*H + 1)); 
SE_CC = SE_CC*pi*E_star/2; 
%
figure; 
hold on
plot(p_bar_dense, SE_Persson01, '-', 'LineWidth', 2);
plot(p_bar_dense, SE_YP08, '-', 'LineWidth', 2);
plot(p_bar_dense, SE_WM17, '-', 'LineWidth', 2);
plot(p_bar_dense, SE_Xu24, '-', 'LineWidth', 2); 
plot(p_bar_dense, SE_Mix, '-', 'LineWidth', 2); 
errorbar(p_bar, Uel_mean, Uel_std, 'o'); 
plot(p_bar_dense, ones(size(p_bar_dense))*SE_CC, '--', 'LineWidth', 2); 
hold off
xlabel('$\bar{p}$ (Pa)', 'interpreter', 'latex'); 
ylabel('Uel (J)');  
legend('Persson01', 'YP08', 'WM17', 'Present work', 'GFMD', 'Present work (hybrid)', 'Uel at complete contact'); 
% save('Fig_R3.mat', 'p_bar_dense', 'SE_Persson01', 'SE_YP08', 'SE_WM17', 'SE_Xu24', 'SE_Mix', ...
%      'p_bar', 'Uel_mean', 'Uel_std', 'SE_CC'); 