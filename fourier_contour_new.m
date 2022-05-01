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
param.z4 = [-0.3, 0.2, 0.4, -0.4]';
param.z5 = [-0.05, 0.08, 0.06, 0.03]';
param.sd = [param.z1; param.z2; param.z3; param.z4; param.z5; 0; 0];


%% simulation
sim.t = 20;
sim.dt = 0.01;
sim.N = 20;

sim.c = zeros(sim.N, 1);
sim.s_step = 1/sim.N;
for i = 1:1:sim.N
    
    s = i*1.5*sim.s_step;
    sim.c(2*i-1) = 8*cos(2*pi*s)*0.8+10;
    sim.c(2*i) = 8*sin(2*pi*s)*0.8+5;
%     sim.c(2*i-1) = normrnd(5, 0.1);
%     sim.c(2*i) = normrnd(5, 0.1);
%     sim.c(2*i-1) = -4.2+0.4*i;
%     sim.c(2*i) = -4.2+0.4*i; 
    
end

sim.step = 2*pi/sim.N;
sim.G = generate_G(sim.step, param.rho-0.5*sim.step, param.N, 1);
sim.k = 2;
sim.L = generate_L(sim.N, sim.k);
sim.L = kron(sim.L, eye(2));

data.u = [];
data.t = [];
data.cd = param.G*param.sd;
data.c = sim.c;
data.ds = [];

ctrl.k = 2;
X0 = sim.c;

for i=0:sim.dt:sim.t
    
    % controller
    [ctrl.u, sim.s, sim.ds] = controller(ctrl.k, sim.L, sim.G, sim.c, param.sd);
    norm(sim.ds)
                                     
    % integration
    tt = [i, i+sim.dt];    
    [~, Temp] = ode45(@(t, X)dynamics(t, X, ctrl.u), tt, X0, ...
        odeset('RelTol', 1e-6, 'AbsTol', 1e-6)); 
    
    % update
    X0 = Temp(end,:)';
    sim.c = X0;
    
    % saving data
    data.t = [data.t, i];
    data.c = [data.c, X0];
    data.u = [data.u, ctrl.u];
    data.ds = [data.ds, sim.ds];
    
end

%% plotting
figure(1)
plot(data.cd(1:2:end), data.cd(2:2:end), 'LineWidth', 1, 'Color', 'k');
hold on; grid on;
plot(data.c(1:2:end, 1), data.c(2:2:end, 1), 'b*',  'LineWidth', 1); hold on;
plot(data.c(1:2:end, end), data.c(2:2:end, end), 'o', 'LineWidth', 1, 'Color', [241, 64, 64]/255); hold on;

for i=1:sim.N

    plot(data.c(2*i-1, :), data.c(2*i, :), 'LineWidth', 1); hold on;
    
end
% set(gca, 'FontSize',13);
xlabel('x positions'); 
ylabel('y positions');
title('Contour-based multi-agent formation');
legend('Desired Contour', 'Agents Initial Position',...
    'Agents Final Position', 'Location', 'NorthWest');

figure(2)
plot(data.cd(1:2:end), data.cd(2:2:end), 'LineWidth', 1, 'Color', 'k');
hold on; grid on;
plot(data.c(1:2:end, end), data.c(2:2:end, end), 'o', 'LineWidth', 1, 'Color', [241, 64, 64]/255 ); hold on;
for i=1:sim.N    
    for j = i-sim.k:i+sim.k
        if j==i
            continue;
        elseif j<=0
            temp = sim.N + j;
        elseif j>sim.N
            temp = j - sim.N;
        else
            temp = j;
        end
        
        plot([data.c(2*i-1,end); data.c(2*temp-1,end)],...
            [data.c(2*i,end); data.c(2*temp,end)],...
            '--', 'LineWidth', 1, 'Color', [26, 111, 223]/255);
        hold on;
    end    
end
% set(gca, 'FontSize',11);
xlabel('x positions'); 
ylabel('y positions');
legend('Desired Contour', 'Agents Final Position',...
    'Links between Neighbors','Location', 'NorthWest');

figure(3)
for i=1:2*sim.N
    
    plot(data.t, data.u(i, :), 'LineWidth', 1); hold on;
    
end
grid on;
xlabel('time(s)'); ylabel('Inputs');

figure(4)
for i=1:size(sim.ds, 1)
    
    plot(data.t, data.ds(i, :), 'LineWidth', 1); hold on;
    
end
grid on;
xlabel('time(s)'); ylabel('ds');

figure(5)
plt.step = 0.01;
plt.G = generate_G(plt.step, param.rho, param.N, 1);
plt.c = plt.G*sim.s;
plot(data.cd(1:2:end), data.cd(2:2:end), 'LineWidth', 1, 'color', 'r');
hold on; 
plot(plt.c(1:2:end), plt.c(2:2:end), '--', 'LineWidth', 1, 'Color', 'k');
grid on;
xlabel('x positions'); ylabel('y positions');
legend('Desired Contour', 'Resulted Contour', 'Location', 'NorthWest');
