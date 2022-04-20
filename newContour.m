clc; clear;
close all;

for t = 0:0.1:200
    
    s = 0:0.001:1;
    
    x = 0.01*(4*t+(2*t*sin(2*pi*s)+800).*cos(2*pi*s));
    y = 0.01*(4*t+(2*t*cos(2*pi*s)+800).*sin(2*pi*s));
    
    plot(x, y, 'k', 'LineWidth', 1);
    grid on;
    axis([-10, 25,-12, 23]);
    drawnow;
    
    t
    
end