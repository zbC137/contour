function [len, len_list] = cLength(start, final, step, t)

s = start:step:final;

x_dot = 0.01*(2*pi*(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800).*sin(2*pi*s)+...
    (-4*pi*t*cos(4*pi*s)+20*pi*t*cos(10*pi*s)).*cos(2*pi*s));
y_dot = 0.01*(-2*pi*(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800).*cos(2*pi*s)+...
    (-4*pi*t*cos(4*pi*s)+20*pi*t*cos(10*pi*s)).*sin(2*pi*s));

len_list = step*sqrt(x_dot.^2+y_dot.^2);

len = sum(len_list);

len_list = len_list/len;

end