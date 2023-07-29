function plotFigure(data, sim, flag, dyn)

id = floor(100/sim.step)+1;
id1 = floor(70/sim.step)+1;
id2 = floor(140/sim.step)+1;

if strcmp(flag, 'fault') || strcmp(flag, 'mixed')
    a = [-10, 25, -10, 25];
else
    a = [-10, 20, -10, 15];
end

figure(2)
set(gcf,'Position', [400, 200, 600, 500], 'Color', 'White');
plot(data.f1.cd(1:2:end,1), data.f1.cd(2:2:end,1), 'k', 'LineWidth', 1.5);
hold on;
plot(data.f1.cd(1:2:end,id), data.f1.cd(2:2:end,id), 'r','LineWidth', 1.5);
hold on;
plot(data.f1.cd(1:2:end,end), data.f1.cd(2:2:end,end), 'b','LineWidth', 1.5);
set(gca, 'FontSize', 12);
legend('Curve at t=0s', 'Curve at t=100s', 'Curve at t=200s',  'Location', 'Best');
% grid on;

figure(3)
set(gcf,'Position', [300, 100, 720, 600], 'Color', 'White');
subplot('Position', [0.07, 0.55 0.4, 0.4])
p1 = plot(data.f1.cd(1:2:end,1), data.f1.cd(2:2:end,1), 'k',...
        data.f1.cds(1:2:end,1), data.f1.cds(2:2:end,1), 'r--',...
        data.f1.cs(1:2:end,1), data.f1.cs(2:2:end,1),'-.',...
        data.c(1:2:end,1), data.c(2:2:end,1), 'b*', 'LineWidth', 1.5);
p1(3).Color = [0.4660 0.95 0.1880];

if strcmp(dyn, 'nonholonomic')
    hold on;
    u = cos(data.theta(:, 1));
    v = sin(data.theta(:, 1));
    quiver(data.c(1:2:end,1), data.c(2:2:end,1), u, v, 1);
end
axis(a);
title('(a) t = 0s')
% grid on;

subplot('Position', [0.55, 0.55, 0.4, 0.4])
p2 = plot(data.f1.cd(1:2:end,id1), data.f1.cd(2:2:end,id1), 'k',...
        data.f1.cds(1:2:end,id1), data.f1.cds(2:2:end,id1), 'r--',...
        data.f1.cs(1:2:end,id1), data.f1.cs(2:2:end,id1),'g-.',...
        data.f1.c(1:2:end, id1), data.f1.c(2:2:end, id1), 'b*', 'LineWidth', 1.5);
p2(3).Color = [0.4660 0.95 0.1880];

if strcmp(dyn, 'nonholonomic')
    hold on;
    u = cos(data.theta(:, id1));
    v = sin(data.theta(:, id1));
    quiver(data.f1.c(1:2:end, id1), data.f1.c(2:2:end, id1), u, v, 0.5);
end
axis(a);
title('(b) t = 70s')
% grid on;

subplot('Position', [0.07, 0.07, 0.4, 0.4])
p3 = plot(data.f1.cd(1:2:end,id2), data.f1.cd(2:2:end,id2), 'k',...
        data.f1.cds(1:2:end,id2), data.f1.cds(2:2:end,id2), 'r--',...
        data.f1.cs(1:2:end,id2), data.f1.cs(2:2:end,id2),'g-.',...
        data.f1.c(1:2:end, id2), data.f1.c(2:2:end, id2), 'b*', 'LineWidth', 1.5);
p3(3).Color = [0.4660 0.95 0.1880];

if strcmp(dyn, 'nonholonomic')
    hold on;
    u = cos(data.theta(:, id2));
    v = sin(data.theta(:, id2));
    quiver(data.f1.c(1:2:end, id2), data.f1.c(2:2:end, id2), u, v, 0.5);
end
axis(a);
title('(c) t = 140s')
% grid on;

subplot('Position', [0.55, 0.07 0.4, 0.4])
% plot(data.f1.cd(1:2:end,end), data.f1.cd(2:2:end,end), 'k',...
%         data.f1.cds(1:2:end,end), data.f1.cds(2:2:end,end), 'r',...
%         data.f1.cs(1:2:end,end), data.f1.cs(2:2:end,end),'g',...
%         data.f1.c(1:2:end, end), data.f1.c(2:2:end, end), 'b*',...
%         data.f1.c(1:2:end, :)', data.f1.c(2:2:end, :)', 'LineWidth', 1);
p4 = plot(data.f1.cd(1:2:end,end), data.f1.cd(2:2:end,end), 'k',...
        data.f1.cds(1:2:end,end), data.f1.cds(2:2:end,end), 'r--',...
        data.f1.cs(1:2:end,end), data.f1.cs(2:2:end,end),'g-.',...
        data.f1.c(1:2:end, end), data.f1.c(2:2:end, end), 'b*', 'LineWidth', 1.5);
p4(3).Color = [0.4660 0.95 0.1880];

if strcmp(dyn, 'nonholonomic')
    hold on;
    u = cos(data.theta(:, end));
    v = sin(data.theta(:, end));
    quiver(data.f1.c(1:2:end, end), data.f1.c(2:2:end, end), u, v, 0.5);
end
axis(a);
title('(d) t = 200s')
% grid on;

figure(4)
set(gcf,'Position', [400, 200, 600, 500], 'Color', 'White');
for i=1:sim.N
    
    plot(data.t, data.u(i, :), 'LineWidth', 1.5, 'Color', 'k'); hold on;
    
end
% grid on;
set(gca, 'FontSize', 12);
xlabel('time(s)'); ylabel('Inputs');

figure(5)
if strcmp(flag, 'switch')
    set(gcf,'Position', [300, 100, 720, 600], 'Color', 'White');
    
    subplot('Position', [0.07, 0.55, 0.4, 0.4])
    plot(graph(-sim.L1+diag(diag(sim.L1))), 'LineWidth', 1.5,  'EdgeColor', 'r',...
        'MarkerSize', 5, 'NodeFontSize', 12);
    title('(a)');
    set(gca, 'FontSize', 12);
    
    subplot('Position', [0.55, 0.55, 0.4, 0.4])
    plot(graph(-sim.L2+diag(diag(sim.L2))), 'LineWidth', 1.5,  'EdgeColor', 'r',...
        'MarkerSize', 5, 'NodeFontSize', 12);
    title('(b)');
    set(gca, 'FontSize', 12);
    
    subplot('Position', [0.3, 0.07, 0.4, 0.4])
    plot(graph(-sim.L3+diag(diag(sim.L3))), 'LineWidth', 1.5,  'EdgeColor', 'r',...
        'MarkerSize', 5, 'NodeFontSize', 12);
    title('(c)'); 
    set(gca, 'FontSize', 12);
else
    set(gcf,'Position', [400, 200, 600, 500], 'Color', 'White');
    plot(graph(-sim.L1+diag(diag(sim.L1))), 'LineWidth', 1.5,  'EdgeColor', 'r',...
        'MarkerSize', 5, 'NodeFontSize', 12);
end

figure(6)
set(gcf,'Position', [400, 200, 600, 500], 'Color', 'White');

subplot('Position', [0.1, 0.55, 0.83, 0.35])
for i=1:size(sim.dcoff, 1)
    
    plot(data.t, data.dcoff(i, :), 'LineWidth', 1.5, 'Color', 'k'); hold on;
    
end
% grid on;
set(gca, 'FontSize', 12);
title('(a) Coefficient Errors');

subplot('Position', [0.1, 0.125, 0.83, 0.3])
plot(data.t, data.e(:, :), 'LineWidth', 1.5, 'Color', 'k');
% grid on;
set(gca, 'FontSize', 12);
xlabel('time(s)'); 
title('(b) Position Errors');

% figure(7)
% set(gcf,'Position', [400, 200, 600, 500], 'Color', 'White');
% plot(data.t, data.xe(:, :), 'LineWidth', 1, 'Color', 'k');
% set(gca, 'FontSize', 15);
% xlabel('time(s)'); ylabel('x_e');
% grid on; 

figure(7)
plot(data.t, data.length(1,:), 'LineWidth', 1, 'Color', 'r');
hold on;
plot(data.t, data.length(2:sim.N+1,:), 'LineWidth', 1, 'Color', 'k');
hold on;
plot(data.t, data.length(sim.N+2:end,:), 'LineWidth', 1, 'Color', 'b');
grid on;

figure(8)
plot(data.t, data.length(1,:), 'LineWidth', 1, 'Color', 'r');
hold on;
plot(data.t, data.length(sim.N+2:end,:), 'LineWidth', 1, 'Color', 'k');
grid on;

if strcmp(dyn, 'nonholonomic')
    figure(9)
    plot(data.t, data.theta(:, 1:size(data.t, 2)), 'LineWidth', 1, 'Color', 'k');
    grid on;
end

end

