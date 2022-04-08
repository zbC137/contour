clc; clear;
close all;

%% Fourier contour
param.N = 5;
param.Ns = 6;
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
param.l = [];
param.Gs = [];

for i=0:0.01:param.rho
    
    l_temp = dist_cal(0, i, param.dist_step, param.N, param.sd(1:20))/param.perimeter;
    param.l = [param.l;l_temp];
    
    Gs_temp = zeros(2, 22);
    for j=1:param.Ns
        
        Gs_temp(:,4*j-3:4*j) = [cos(2*pi*j*l_temp), sin(2*pi*j*l_temp),        0,       0
                                           0,        0,         cos(2*pi*j*l_temp), sin(2*pi*j*l_temp)];
                                    
    end
    
    Gs_temp(:, 4*param.Ns+1:4*param.Ns+2) = eye(2);
    param.Gs = [param.Gs; Gs_temp];
    
    i
    
end


data.cd = param.G*param.sd;

%% curve fitting-Newton method
param.coff = zeros(4*param.Ns+2,1)+1;
% param.coff = pinv(param.Gs'*param.Gs)*param.Gs'*data.cd;
ndelta = 1e+4;

index = 1;
fs = param.Gs*param.coff;
e = fs-data.cd;
while ndelta>1e-4
    
    index    
    
    delta = -pinv(param.Gs'*param.Gs)*param.Gs'*e;
    param.coff = param.coff+delta;
    fs = param.Gs*param.coff;
    e = fs-data.cd;
    error = norm(e);
    ndelta = norm(delta)
    
    index = index+1;
    
end

data.cds = param.Gs*param.coff;

%% formation control simulation
sim.t = 5;
sim.dt = 0.001;
sim.N = 10;

sim.c = zeros(sim.N, 1);
for i = 1:2*sim.N
    
    if mod(i, 2)==0
        sim.c(i) = normrnd(param.coff(4*param.Ns+2),0.1);
    else
        sim.c(i) = normrnd(param.coff(4*param.Ns+1),0.1);
    end
    
end

sim.step = 1/sim.N;
sim.Gs = generate_Gs(sim.step, 1-0.5*sim.step, param.Ns);
sim.k = 2;
if sim.k == 0
    sim.L = eye(2*sim.N);
else
    sim.L = eye(2*sim.N)*sim.k*2;
end

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
        
        sim.L(2*i-1,2*temp-1) = -1;
        sim.L(2*i, 2*temp) = -1;
    end
end

data.ref = sim.Gs*param.coff;
data.u = [];
data.t = [];
data.c = sim.c;
data.dcoff = [];

ctrl.k = 2;
X0 = sim.c;

for i=0:sim.dt:sim.t
    
    % controller
    [ctrl.u, sim.coff, sim.dcoff] = controller(ctrl.k, sim.L, sim.Gs,...
                                         4*param.Ns+2, sim.c, param.coff);
    norm(sim.dcoff);
                                     
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
    data.dcoff = [data.dcoff, sim.dcoff];
    
    i
    
end

%% plotting
figure(1)
plot(data.cd(1:2:end), data.cd(2:2:end), 'LineWidth', 1, 'Color', 'k');
hold on; grid on;
plot(data.cds(1:2:end), data.cds(2:2:end), 'LineWidth', 1, 'Color', 'r');
plot(data.c(1:2:end, 1), data.c(2:2:end, 1), 'b*',  'LineWidth', 1); hold on;
plot(data.ref(1:2:end), data.ref(2:2:end), 'ko', 'LineWidth', 1); hold on;
plot(data.c(1:2:end, end), data.c(2:2:end, end), 'o', 'LineWidth', 1, 'Color', [241, 64, 64]/255); hold on;

for i=1:sim.N

    plot(data.c(2*i-1, :), data.c(2*i, :), 'LineWidth', 1); hold on;
    
end
% set(gca, 'FontSize',13);
xlabel('x positions'); 
ylabel('y positions');
% title('Contour-based multi-agent formation');
legend('Desired Contour', 'Re-parameterized Contour', 'Agents Initial Position',...
    'Agents Target Position', 'Agents Final Position', 'Location', 'NorthWest');

figure(2)
plot(data.cds(1:2:end), data.cds(2:2:end), 'LineWidth', 1, 'Color', 'k');
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
legend('Re-parameterized Contour', 'Agents Final Position',...
    'Links between Neighbors','Location', 'NorthWest');

figure(3)
for i=1:2*sim.N
    
    plot(data.t, data.u(i, :), 'LineWidth', 1); hold on;
    
end
grid on;
xlabel('time(s)'); ylabel('Inputs');

figure(4)
for i=1:size(sim.dcoff, 1)
    
    plot(data.t, data.dcoff(i, :), 'LineWidth', 1); hold on;
    
end
grid on;
xlabel('time(s)'); ylabel('coefficient errors');

figure(5)
plt.step = 0.001;
plt.Gs = generate_Gs(plt.step, 1, param.Ns);
plt.c = plt.Gs*sim.coff;
plot(data.cds(1:2:end), data.cds(2:2:end), 'LineWidth', 1, 'color', 'r');
hold on; 
plot(plt.c(1:2:end), plt.c(2:2:end), '--', 'LineWidth', 1, 'Color', 'k');
grid on;
xlabel('x positions'); ylabel('y positions');
legend('Re-parameterized Contour', 'Resulted Contour', 'Location', 'NorthWest');
