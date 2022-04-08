clc; clear;
close all;

%% parameters
param.N = 5;
param.rho = 2*pi;
param.step = 0.01;
param.dist_step = 0.001;
param.G = generate_G(param.step, param.rho, param.N, 1);

param.z1 = [2, -4, 4, 2]';
param.z2 = [-1.3, 1.2, 1.4, -1.4]';
param.z3 = [0.8, -0.7, 0.9, -0.9]';
param.z4 = [-0.3, 0.2, 0.4, -0.4]';
param.z5 = [-0.05, 0.08, 0.06, 0.03]';
param.sd = [param.z1; param.z2; param.z3; param.z4; param.z5; 0; 0];
param.perimeter = dist_cal(0, param.rho, param.dist_step,...
    param.N, param.sd(1:20));

%% simulation
sim.t = 20;
sim.dt = 0.01;
sim.N = 12;

sim.c = zeros(sim.N, 1);
for i = 1:2*sim.N
    
    sim.c(i) = normrnd(0,0.1);
    
end

sim.step = 2*pi/sim.N;
sim.rho = [];
for i=1:sim.N
    
    sim.rho = [sim.rho; sim.step*(i-1)];
    
end
% sim.G = generate_G(sim.step, param.rho-0.5*sim.step, param.N);
sim.G = generate_G(sim.N, sim.rho, param.N, 2);

sim.dist_ref = zeros(sim.N,1)+param.perimeter/sim.N;
sim.k = 1;
if sim.k == 0
    sim.L1 = eye(2*sim.N);
else
    sim.L1 = eye(2*sim.N)*sim.k*2;
end

sim.L2 = eye(sim.N)*sim.k*2;

for i = 1:sim.N
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
        
        sim.L1(2*i-1,2*temp-1) = -1;
        sim.L1(2*i, 2*temp) = -1;
        
        sim.L2(i, temp) = -1;        
        
    end
end

sim.P = eye(sim.N) - sim.dt*sim.L2;
sim.dists = zeros(sim.N, 1);

data.u = [];
data.t = [];
data.cd = param.G*param.sd;
data.c = sim.c;
data.ds = [];

ctrl.k = 2;
ctrl.p = 0.1;
X0 = sim.c;

for i=0:sim.dt:sim.t
    i
    
    % controller
    [ctrl.u, sim.s, sim.ds] = controller(ctrl.k, sim.L1, sim.G,...
                                         4*param.N+2, sim.c, param.sd);
    
    % integration
    tt = [i, i+sim.dt];    
    [~, Temp] = ode45(@(t, X)dynamics(t, X, ctrl.u), tt, X0, ...
        odeset('RelTol', 1e-6, 'AbsTol', 1e-6));
    
    % distances
    for j=1:sim.N-1        
        sim.dists(j) = dist_cal(sim.rho(j), sim.rho(j+1), ...
            param.dist_step, param.N, sim.s(1:20));        
    end
    sim.dists(end) = dist_cal(sim.rho(end), sim.rho(1)+2*pi,...
        param.dist_step, param.N, sim.s(1:20));
    
    % update
    X0 = Temp(end,:)';
    sim.c = X0;
%     sim.rho = sim.rho - ctrl.p*sim.P*sim.dists - ctrl.p*sim.dt*sim.L2*sim.dist_ref;
    sim.rho = sim.rho + ctrl.p*(sim.dists-[sim.dists(end); sim.dists(1:end-1)]);
    sim.G = generate_G(sim.N, sim.rho, param.N, 2);
    
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
plot(data.c(1:2:end, 1), data.c(2:2:end, 1), 'b*'); hold on;
plot(data.c(1:2:end, end), data.c(2:2:end, end), 'ro'); hold on;

for i=1:sim.N
    
    plot(data.c(2*i-1, :), data.c(2*i, :), 'LineWidth', 1); hold on;
    
end

xlabel('x positions'); ylabel('y positions');
title('Contour-based multi-agent formation');
legend('Desired Contour', 'Agents Initial Position',...
    'Agents Final Position', 'Location', 'NorthWest');

figure(2)
for i=1:2*sim.N
    
    plot(data.t, data.u(i, :), 'LineWidth', 1); hold on;
    
end
grid on;
xlabel('time(s)'); ylabel('Inputs');

figure(3)
for i=1:size(sim.ds, 1)
    
    plot(data.t, data.ds(i, :), 'LineWidth', 1); hold on;
    
end
grid on;
xlabel('time(s)'); ylabel('ds');

figure(4)
plt.step = 0.01;
plt.G = generate_G(plt.step, param.rho, param.N, 1);
plt.c = plt.G*sim.s;
plot(data.cd(1:2:end), data.cd(2:2:end), 'LineWidth', 1, 'color', 'r');
hold on; 
plot(plt.c(1:2:end), plt.c(2:2:end), '--', 'LineWidth', 1, 'Color', 'k');
grid on;
xlabel('x positions'); ylabel('y positions');
legend('Desired Contour', 'Resulted Contour', 'Location', 'NorthWest');
