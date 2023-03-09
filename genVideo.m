function genVideo(dat, data)
% clear; clc;

%% Data processing
N = size(data.t, 2);
data.dcoff = abs(data.dcoff);
data.e = abs(data.e);
data.theta = data.theta(:, 1:N);

%% Parameters
rate = 30;
quality = 100;

if strcmp(dat, 'fault')
    a = [-10, 25, -10, 25];
else
    a = [-10, 20, -10, 15];
end

figure()

aviobj=VideoWriter(['video\', dat, '.mp4'], 'MPEG-4');%新建叫example.avi的文件
aviobj.FrameRate = rate;           % frame rate 帧率
aviobj.Quality = quality;
open(aviobj); %    打开对象

for i = 1:5:N
    set(gcf, 'Color', 'White');
    plot(data.f1.cd(1:2:end,i), data.f1.cd(2:2:end,i), 'k',...
        data.f1.cds(1:2:end,i), data.f1.cds(2:2:end,i), 'r--',...
        data.c(1:2:end,i), data.c(2:2:end,i), 'b*', 'LineWidth', 1.5);
    hold on;
    u = cos(data.theta(:, i));
    v = sin(data.theta(:, i));
    quiver(data.c(1:2:end,i), data.c(2:2:end,i), u, v, 0.2, ...,
        'Color', '#77AC30', 'LineWidth', 1.5);
    axis(a); 
    set(gca, 'FontSize', 12);
    title('Simulation Results');
    legend('Observed Curve', 'Approximated Curve', 'Agents', 'Orientations',...
        'Location', 'best', 'NumColumns', 2, 'FontSize', 9);
    
    currFrame = getframe(gcf);   % 获取当前帧
    writeVideo(aviobj,currFrame);   % 保存当前帧
    delete(gca);
    i
end
close(aviobj);   % 关闭保存视频

aviobj=VideoWriter(['video\', dat, 'Coff.mp4'], 'MPEG-4');%新建叫example.avi的文件
aviobj.FrameRate = rate;           % frame rate 帧率
aviobj.Quality = quality;
open(aviobj); %    打开对象
for i=1:5:N        % 生成每帧图像 一共N帧
    set(gcf, 'Color', 'White');
    if strcmp(dat, 'less')
        plot(data.t(1:i), data.lim(1:i), 'LineWidth', 1.5, 'Color', 'b');
        set(gca, 'FontSize', 12);
        title('Coefficient Error Metric')
    else
        plot(data.t(1:i), data.dcoff(:, 1:i), 'LineWidth', 1.5, 'Color', 'b');
        set(gca, 'FontSize', 12);
        title('Coefficient Errors')
    end
    xlabel('time (s)')
   
    currFrame = getframe(gcf);   % 获取当前帧
    writeVideo(aviobj,currFrame);   % 保存当前帧
    delete(gca);
    i
end
close(aviobj); 

aviobj=VideoWriter(['video\', dat, 'Pos.mp4'], 'MPEG-4');%新建叫example.avi的文件
aviobj.FrameRate = rate;           % frame rate 帧率
aviobj.Quality = quality;
open(aviobj); %    打开对象
for i=1:5:N        % 生成每帧图像 一共N帧
    set(gcf, 'Color', 'White');      
    plot(data.t(1:i), data.e(1, 1:i), 'LineWidth', 1.5, 'Color', 'r'); 
    hold on;
    plot(data.t(1:i), data.e(2, 1:i), 'LineWidth', 1.5, 'Color', '#77AC30'); 
    hold on;
    plot(data.t(1:i), data.e(3:2:end, 1:i), 'LineWidth', 1.5, 'Color', 'r'); 
    hold on;
    plot(data.t(1:i), data.e(4:2:end, 1:i), 'LineWidth', 1.5, 'Color', '#77AC30');
    set(gca, 'FontSize', 12);
    xlabel('time (s)');
    title('Position Errors');
    legend('X-axis', 'Y-axis', 'Location', 'best','NumColumns', 2, 'FontSize', 9)
   
    currFrame = getframe(gcf);   % 获取当前帧
    writeVideo(aviobj,currFrame);   % 保存当前帧
    delete(gca);
    i
end
close(aviobj); 

aviobj=VideoWriter(['video\', dat, 'Ori.mp4'], 'MPEG-4');%新建叫example.avi的文件
aviobj.FrameRate = rate;           % frame rate 帧率
aviobj.Quality = quality;
open(aviobj); %    打开对象
for i=1:5:N        % 生成每帧图像 一共N帧
    set(gcf, 'Color', 'White');      
    plot(data.t(1:i), data.theta(:, 1:i), 'LineWidth', 1.5, 'Color', 'k'); 
    set(gca, 'FontSize', 12);
    xlabel('time (s)');
    title('Orientations');
   
    currFrame = getframe(gcf);   % 获取当前帧
    writeVideo(aviobj,currFrame);   % 保存当前帧
    delete(gca);
    i
end
close(aviobj); 

end