clc; clear;
close all;

%% parameters
param.N = 5;
param.rho = 2*pi;
param.step = 0.01;
param.G = generate_G(param.step, param.rho, param.N, 1);

param.z1 = [2, -4, 4, 2]';
param.z2 = [-1.3, 1.2, 1.4, -1.4]';
param.z3 = [0.8, -0.7, 0.9, -0.9]';
param.z4 = [-0.3,  0.2, 0.4, -0.4]';
param.z5 = [-0.05, -0.08, -0.06, -0.03]';
param.sd = [param.z1; param.z2; param.z3; param.z4; param.z5; 0; 0];

%% plotting
data.c = param.G*param.sd;

for i=0.1:0.01:1.5
    
    plot(i*data.c(1:2:end), i*data.c(2:2:end), 'LineWidth', 1, 'Color', 'k');
    axis([-11, 9, -8, 12]);
    pause(0.05);
    drawnow;

end
