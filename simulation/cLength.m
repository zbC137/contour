function [len, len_list] = cLength(start, final, step, t, option)

s = start:step:final;

if strcmp(option, 'fault')
    x = 0.01*(4*t+(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800).*cos(2*pi*s));
    y = 0.01*(2*t+(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800).*sin(2*pi*s));
    
    % x_dot = 0.01*(-2*pi*(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800).*sin(2*pi*s)+...
    %     (4*pi*t*cos(4*pi*s)-20*pi*t*cos(10*pi*s)).*cos(2*pi*s));
    % y_dot = 0.01*(2*pi*(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800).*cos(2*pi*s)+...
    %     (4*pi*t*cos(4*pi*s)-20*pi*t*cos(10*pi*s)).*sin(2*pi*s));
elseif strcmp(option, 'mixed')
    if t<100
%     if t<6
        x = 0.01*(4*t+(t*sin(4*pi*s)+800)*cos(2*pi*s));
        y = 0.01*(4*t+(t*cos(4*pi*s)+800)*sin(2*pi*s));
%         x = -7+10*cos(2*pi*s);
%         y = 10*sin(2*pi*s);
%     elseif t<12
%         x = 8*cos(2*pi*s);
%         y = 12*sin(2*pi*s);
    else
        x = 0.01*(4*(t-100)+((t-100)*sin(4*pi*s)+2*(t-100)*cos(10*pi*s)+800)*cos(2*pi*s));
        y = 0.01*(2*(t-100)+((t-100)*sin(4*pi*s)+2*(t-100)*cos(10*pi*s)+800)*sin(2*pi*s));
%         x = 4+(8*(1-sin(2*pi*s))).*cos(2*pi*s);
%         y = 4+ (8*(1-sin(2*pi*s))).*sin(2*pi*s);
    end
else
    x = 0.01*(4*t+(t*sin(4*pi*s)+800).*cos(2*pi*s));
    y = 0.01*(4*t+(t*cos(4*pi*s)+800).*sin(2*pi*s));
%     x = 2*(2-2.^(sin(5*2*pi*s))).*cos(2*pi*s);
%             y = 2*(2-2.^(sin(5*2*pi*s))).*sin(2*pi*s);
% x = 15*sin(t*cos(t)*2*pi*s);
%                y = 15*sin(t*sin(t)*2*pi*s);
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