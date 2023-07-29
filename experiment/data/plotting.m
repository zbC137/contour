clear; clc;
close all;

load('data4.mat');
an = 9;
cn = 6;

t = data(1:end-1,2*an+1);
coff_err = data(2:end, 2*an+2:2*an+4*cn+3)';
pos_err = data(2:end,2*an+4*cn+4:4*an+4*cn+3)';
ori = data(2:end, 4*an+4*cn+4:end);

dcoff = [t, vecnorm(coff_err)'*1.5e-3];
e = [t, vecnorm(pos_err)'*1.5e-3];
theta = [t, ori];
xi = [t, data(2:end, end)*1.5e-3];

figure(1)
plot(t, pos_err*1.5e-3, 'LineWidth', 1, 'Color', 'k');
title("Position Errors");
xlabel("time (s)"); ylabel("pixels");

figure(2)
plot(t, coff_err*1.5e-3, 'LineWidth', 1, 'Color', 'k');
title("Coefficient Errors");
xlabel("time (s)");

figure(3)
plot(t, ori, 'LineWidth', 1, 'Color', 'k');
title("Orientations");
xlabel("time (s)"); ylabel("rad");

