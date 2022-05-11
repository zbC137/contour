clc; clear;
close all;

tic;
%% parameters
param.Ns = 6;
% param.Ns = 3;
param.step= 0.001;
% param.Gs = generate_Gs(param.step, 1, param.Ns);

%% simulation
sim.N = 8;     % agent number
sim.t = 200;    % simulation time
sim.step = 0.1;
sim.s_step = 1/sim.N;

sim.c = zeros(2*sim.N, 1);      % initial positions
for i = 1:1:sim.N
    
%     s = i*sim.s_step;
%     sim.c(2*i-1) = 8*cos(2*pi*s)*0.2+10;
%     sim.c(2*i) = 8*sin(2*pi*s)*0.2;
    sim.c(2*i-1) = normrnd(5, 2);
    sim.c(2*i) = normrnd(5, 2);
%     sim.c(2*i-1) = -4.2+0.4*i;
%     sim.c(2*i) = -4.2+0.4*i; 
    
end
sim.c0 = sim.c;

sim.Gs = generate_Gs(0, sim.s_step, 1-sim.s_step/2, param.Ns);

% Laplacian
sim.L2 = generate_L(sim.N, 2);
sim.L = kron(sim.L2, eye(2));

% data recording
% data.ref = sim.Gs*param.coff;
data.u = [];
data.t = [];
data.c = sim.c;
data.dcoff = [];
data.f1.cd = [];
data.f1.cds = [];
data.f1.cs = [];
data.f1.c = [];
data.c0 = sim.c;
data.e = [];
data.xe = [];

% controller parameters
ctrl.k = 2;
X0 = sim.c;

% show the real-time simulation
figure(1)

for t = 0:sim.step:sim.t
       
    data.cd = [];
    
    for s = 0:param.step:1
        
%         x = 0.01*(4*t+(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800)*cos(2*pi*s));
%         y = 0.01*(2*t+(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800)*sin(2*pi*s));

        x = 0.01*(4*t+(t*sin(4*pi*s)+800)*cos(2*pi*s));
        y = 0.01*(4*t+(t*cos(4*pi*s)+800)*sin(2*pi*s));
        
        data.cd = [data.cd; [x; y]];
        
    end
    
    [len, len_list] = cLength(0, 1, param.step, t);
    param.Gs = generate_Gs_new(len_list, param.Ns);
    
    % curve fitting
    param.coff = zeros(4*param.Ns+2,1)+1;
    % param.coff = pinv(param.Gs'*param.Gs)*param.Gs'*data.cd;
    ndelta = 1e+4;
    
    index = 1;
    fs = param.Gs*param.coff;
    e = fs-data.cd;
    while ndelta>1e-4
        
        index;
        
        delta = -pinv(param.Gs'*param.Gs)*param.Gs'*e;
        param.coff = param.coff+delta;
        fs = param.Gs*param.coff;
        e = fs-data.cd;
        error = norm(e);
        ndelta = norm(delta);
        
        index = index+1;
        
    end
    data.cds = param.Gs*param.coff;
    
    data.e = [data.e, sim.c-sim.Gs*param.coff];
    data.xe = [data.xe, sim.Gs*pinv(sim.Gs)*sim.c-sim.Gs*param.coff];
    
%     sim.Gs = generate_Gs(0.005*t, sim.s_step, 0.005*t+1-sim.s_step/2, param.Ns);

    % control
     [ctrl.u, sim.coff, sim.dcoff] = controller(ctrl.k, sim.L, sim.Gs, sim.c, param.coff);
        
    % integration
    tt = [t, t+sim.step];    
    [~, Temp] = ode45(@(t, X)dynamics(t, X, ctrl.u), tt, X0, ...
        odeset('RelTol', 1e-6, 'AbsTol', 1e-6));
    
    % update
    X0 = Temp(end,:)';
    sim.c = X0;
    
    % saving data
    data.t = [data.t, t];
    data.c = [data.c, X0];
    data.u = [data.u, ctrl.u];
    data.dcoff = [data.dcoff, sim.dcoff];
    data.cs = param.Gs*sim.coff;
    
    data.f1.cd = [data.f1.cd, data.cd];
    data.f1.cds = [data.f1.cds, data.cds];
    data.f1.c = [data.f1.c, data.c(:,end)];
    data.f1.cs = [data.f1.cs, data.cs];
    
    % plot real-time data        
    plot(data.cd(1:2:end), data.cd(2:2:end), 'k',...
        data.cds(1:2:end), data.cds(2:2:end), 'r',...
        data.cs(1:2:end), data.cs(2:2:end),'g',...
        data.c(1:2:end, end), data.c(2:2:end, end), 'b*', 'LineWidth', 1);
    grid on;
    axis([-10, 25,-12, 23]);
    legend('Original Curve', 'Re-paramerterized Curve', 'Real-time Curve',...
        'Agents Real-time Positions',  'Location', 'NorthWest');
    drawnow;
    
    t
    
end

toc;

%% plotting
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
axis([-10,20,-10,15]);
grid on;

figure(32)
plot(data.f1.cd(1:2:end,701), data.f1.cd(2:2:end,701), 'k',...
        data.f1.cds(1:2:end,701), data.f1.cds(2:2:end,701), 'r',...
        data.f1.cs(1:2:end,701), data.f1.cs(2:2:end,701),'g',...
        data.f1.c(1:2:end, 701), data.f1.c(2:2:end, 701), 'b*', 'LineWidth', 1);
xlabel('x positions'); ylabel('y positions');
axis([-10,20,-10,15]);
grid on;

figure(33)
plot(data.f1.cd(1:2:end,1401), data.f1.cd(2:2:end,1401), 'k',...
        data.f1.cds(1:2:end,1401), data.f1.cds(2:2:end,1401), 'r',...
        data.f1.cs(1:2:end,1401), data.f1.cs(2:2:end,1401),'g',...
        data.f1.c(1:2:end, 1401), data.f1.c(2:2:end, 1401), 'b*', 'LineWidth', 1);
xlabel('x positions'); ylabel('y positions');
axis([-10,20,-10,15]);
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
axis([-10,20,-10,15]);
grid on;

% axis([-10, 25,-12, 23]);
% xlabel('x positions'); ylabel('y positions');
% legend('Original Curve', 'Re-paramerterized Curve', 'Real-time Curve',...
%     'Agents Real-time Positions',  'Location', 'NorthWest');

% figure(4)
% plot(data.cds(1:2:end), data.cds(2:2:end), 'LineWidth', 1, 'Color', 'k');
% hold on; grid on;
% plot(data.c(1:2:end, end), data.c(2:2:end, end), 'o', 'LineWidth', 1, 'Color', [241, 64, 64]/255 ); hold on;
% for i=1:sim.N    
%     for j = i-sim.k:i+sim.k
%         if j==i
%             continue;
%         elseif j<=0
%             temp = sim.N + j;
%         elseif j>sim.N
%             temp = j - sim.N;
%         else
%             temp = j;
%         end
%         
%         plot([data.c(2*i-1,end); data.c(2*temp-1,end)],...
%             [data.c(2*i,end); data.c(2*temp,end)],...
%             '--', 'LineWidth', 1, 'Color', [26, 111, 223]/255);
%         hold on;
%     end    
% end
% % set(gca, 'FontSize',11);
% xlabel('x positions'); 
% ylabel('y positions');
% legend('Re-parameterized Contour', 'Agents Final Position',...
%     'Links between Neighbors','Location', 'NorthWest');

figure(5)
for i=1:2*sim.N
    
    plot(data.t, data.u(i, :), 'LineWidth', 1); hold on;
    
end
grid on;
xlabel('time(s)'); ylabel('Inputs');

figure(6)
for i=1:size(sim.dcoff, 1)
    
    plot(data.t, data.dcoff(i, :), 'LineWidth', 1); hold on;
    
end
grid on;
xlabel('time(s)'); ylabel('coefficient errors');

figure(7)
plot(graph(-sim.L2+diag(diag(sim.L2))));

figure(8)
plot(data.t, data.e(1:2:end, :), 'LineWidth', 1);
grid on;
xlabel('time(s)'); ylabel('x-position errors');

figure(9)
plot(data.t, data.e(2:2:end, :), 'LineWidth', 1);
xlabel('time(s)'); ylabel('y-position errors');
grid on;

figure(10)
plot(data.t, data.xe(1:2:end, :), 'LineWidth', 1);
grid on;
xlabel('time(s)'); ylabel('x-position errors fake');

figure(11)
plot(data.t, data.xe(2:2:end, :), 'LineWidth', 1);
xlabel('time(s)'); ylabel('y-position errors fake');
grid on; 
