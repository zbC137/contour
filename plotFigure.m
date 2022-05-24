function plotFigure(data, sim, flag)

if strcmp(flag, 'fault')
    a = [-10, 25, -10, 25];
else
    a = [-10, 20, -10, 15];

figure(2)
plot(data.f1.cd(1:2:end,1), data.f1.cd(2:2:end,1), 'k', 'LineWidth', 1);
hold on;
plot(data.f1.cd(1:2:end,1001), data.f1.cd(2:2:end,1001), 'r','LineWidth', 1);
hold on;
plot(data.f1.cd(1:2:end,end), data.f1.cd(2:2:end,end), 'b','LineWidth', 1);
xlabel('x positions'); ylabel('y positions');
legend('Curve at t=0s', 'Curve at t=100s', 'Curve at t=200s',  'Location', 'Best');
grid on;

figure(31)
plot(data.f1.cd(1:2:end,1), data.f1.cd(2:2:end,1), 'k',...
        data.f1.cds(1:2:end,1), data.f1.cds(2:2:end,1), 'r',...
        data.f1.cs(1:2:end,1), data.f1.cs(2:2:end,1),'g',...
        data.c0(1:2:end), data.c0(2:2:end), 'b*', 'LineWidth', 1);
xlabel('x positions'); ylabel('y positions');
axis(a);
grid on;

figure(32)
plot(data.f1.cd(1:2:end,701), data.f1.cd(2:2:end,701), 'k',...
        data.f1.cds(1:2:end,701), data.f1.cds(2:2:end,701), 'r',...
        data.f1.cs(1:2:end,701), data.f1.cs(2:2:end,701),'g',...
        data.f1.c(1:2:end, 701), data.f1.c(2:2:end, 701), 'b*', 'LineWidth', 1);
xlabel('x positions'); ylabel('y positions');
axis(a);
grid on;

figure(33)
plot(data.f1.cd(1:2:end,1401), data.f1.cd(2:2:end,1401), 'k',...
        data.f1.cds(1:2:end,1401), data.f1.cds(2:2:end,1401), 'r',...
        data.f1.cs(1:2:end,1401), data.f1.cs(2:2:end,1401),'g',...
        data.f1.c(1:2:end, 1401), data.f1.c(2:2:end, 1401), 'b*', 'LineWidth', 1);
xlabel('x positions'); ylabel('y positions');
axis(a);
grid on;

figure(34)
% plot(data.f1.cd(1:2:end,end), data.f1.cd(2:2:end,end), 'k',...
%         data.f1.cds(1:2:end,end), data.f1.cds(2:2:end,end), 'r',...
%         data.f1.cs(1:2:end,end), data.f1.cs(2:2:end,end),'g',...
%         data.f1.c(1:2:end, end), data.f1.c(2:2:end, end), 'b*',...
%         data.f1.c(1:2:end, :)', data.f1.c(2:2:end, :)', 'LineWidth', 1);
plot(data.f1.cd(1:2:end,end), data.f1.cd(2:2:end,end), 'k',...
        data.f1.cds(1:2:end,end), data.f1.cds(2:2:end,end), 'r',...
        data.f1.cs(1:2:end,end), data.f1.cs(2:2:end,end),'g',...
        data.f1.c(1:2:end, end), data.f1.c(2:2:end, end), 'b*', 'LineWidth', 1);
xlabel('x positions'); ylabel('y positions');
axis(a);
grid on;

figure(4)
for i=1:2*sim.N
    
    plot(data.t, data.u(i, :), 'LineWidth', 1, 'Color', 'k'); hold on;
    
end
grid on;
xlabel('time(s)'); ylabel('Inputs');

figure(5)
for i=1:size(sim.dcoff, 1)
    
    plot(data.t, data.dcoff(i, :), 'LineWidth', 1, 'Color', 'k'); hold on;
    
end
grid on;
xlabel('time(s)'); ylabel('coefficient errors', 'Color', 'k');

figure(6)
if strcmp(flag, 'switch')
    subplot(2,2,1)
    plot(graph(-sim.L1+diag(diag(sim.L1))));
    title('(a)');
    subplot(2,2,2)
    plot(graph(-sim.L2+diag(diag(sim.L2))));
    title('(b)');
    subplot(2,2,3)
    plot(graph(-sim.L3+diag(diag(sim.L3))));
    title('(c)');
else
    plot(graph(-sim.L1+diag(diag(sim.L1))));
end

figure(7)
plot(data.t, data.e(:, :), 'LineWidth', 1, 'Color', 'k');
grid on;
xlabel('time(s)'); ylabel('position errors');

figure(8)
plot(data.t, data.xe(:, :), 'LineWidth', 1, 'Color', 'k');
xlabel('time(s)'); ylabel('x_e');
grid on; 

end

