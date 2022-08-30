% Analyze transfer functions, pole-zero plots, and frequency responses of
% all system iterations (with and without hand impedance, with and without
% feedforward force).

clear all;
close all;
clc;

% Define variables
% m = 1.1; % Mass of ball (kg)
% M = 1.9; % Mass of cup (kg)
% l = 0.5; % Length of the massless pendulum rod (m)
% g = 9.8; % Gravity (m/s^2)
% K = 100; % N/m
% B = 10;  % N-s/m

m = 1; % Mass of ball (kg)
M = 1; % Mass of cup (kg)
l = 0.5; % Length of the massless pendulum rod (m)
g = 10; % Gravity (m/s^2)
K = 100; % N/m
B = 10;  % N-s/m

%% No hand impedance
% Transfer function from x (cart position) to theta (ball angle)
% x_to_theta = tf([-1/l 0 0], [1 0 g/l]);
% 
% figure()
% bode(x_to_theta)
% title("Transfer Function: x to {\theta}")
% 
% figure()
% pzmap(x_to_theta)
% title("Transfer Function: x to {\theta}")

%%
% Transfer function from x_0 (ZFT) to x (cart position)

xo_to_x = tf([B K B*g/l K*g/l], [M B ((m+M)*g/l+K) B*g/l K*g/l]);

% figure()
% bode(xo_to_x)
% title("Transfer Function: x_o to x")

figure()
pzmap(xo_to_x)
title("Transfer Function: x_o to x")