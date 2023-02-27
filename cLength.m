function [len, len_list] = cLength(start, final, step, t, option)

s = start:step:final;

if strcmp(option, 'fault')
    x = 0.01*(4*t+(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800).*cos(2*pi*s));
    y = 0.01*(2*t+(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800).*sin(2*pi*s));
    
    % x_dot = 0.01*(-2*pi*(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800).*sin(2*pi*s)+...
    %     (4*pi*t*cos(4*pi*s)-20*pi*t*cos(10*pi*s)).*cos(2*pi*s));
    % y_dot = 0.01*(2*pi*(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800).*cos(2*pi*s)+...
    %     (4*pi*t*cos(4*pi*s)-20*pi*t*cos(10*pi*s)).*sin(2*pi*s));
    
else
    x = 0.001*(4*t+(t*sin(4*pi*s)+800).*cos(2*pi*s));
    y = 0.001*(4*t+(t*cos(4*pi*s)+800).*sin(2*pi*s));
    
    % x_dot = 0.01*(-2*pi*(t*sin(4*pi*s)+800).*sin(2*pi*s)+...
    %     4*pi*t*cos(3*pi*s).*cos(2*pi*s));
    % y_dot = 0.01*(2*pi*(t*cos(4*pi*s)+800).*cos(2*pi*s)-...
    %     4*pi*t*sin(4*pi*s).*sin(2*pi*s));
end

% len_list = step*sqrt(x_dot.^2+y_dot.^2);
len_list = sqrt((x(2:end)-x(1:1000)).^2+(y(2:end)-y(1:1000)).^2);
len_list = [0, len_list];

len = sum(len_list);

len_list = len_list/len;

end