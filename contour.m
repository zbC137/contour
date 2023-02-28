function [data, sim] = contour(option, dys, Ns, N)

tic;
%% parameters
param.Ns = Ns;
% param.Ns = 3;
param.step= 0.001;
% param.Gs = generate_Gs(param.step, 1, param.Ns);

%% simulation
sim.N = N;     % agent number
sim.t = 200;    % simulation time
sim.step = 0.1;
sim.s_step = 1/sim.N;

sim.c = zeros(2*sim.N, 1);      % initial positions and orientations
for i = 1:1:sim.N
    
%         s = i*sim.s_step;
%         sim.c(2*i-1) = (8*cos(2*pi*s)*1.1+10)*0.1;
%         sim.c(2*i) = 8*sin(2*pi*s)*1.1*0.1;
    %     sim.c(2*i-1) = normrnd(5, 2);
    %     sim.c(2*i) = normrnd(5, 2);
    sim.c(2*i-1) = 13;
    sim.c(2*i) = (-8.4+0.5*i);
    
end

sim.Gs = generate_Gs(0, sim.s_step, 1-sim.s_step/2, param.Ns);

% Laplacian
if strcmp(option, 'switch')
    sim.L1 = generate_L(sim.N, 2);
    sim.L2= [2, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        -1, 3, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        -1, -1, 4, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        0, -1, -1, 4, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0
        0, 0, -1, -1, 4, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0
        0, 0, 0, -1, -1, 3, -1, 0, 0, 0, 0, 0, 0, 0, 0
        0, 0, 0, 0, -1, -1, 2, 0, 0, 0, 0, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];
    sim.L3 = generate_L(sim.N, 0);
else
    sim.L1 = generate_L(sim.N, 2);
end

% curve length recording
param.G_cl = generate_Gs(0, 1/2400, 1, param.Ns);
% size(param.G_cl);

if strcmp(dys, 'nonholonomic')
    param.l = 0.1;
    sim.theta = zeros(sim.N, 1) + pi;
    data.theta = sim.theta;
    X0 = [sim.c; sim.theta];
    for i=1:sim.N
        sim.c(2*i-1) = sim.c(2*i-1)+0.01*cos(sim.theta(i));
        sim.c(2*i) = sim.c(2*i)+0.01*sin(sim.theta(i));
    end
else
    X0 = sim.c;
end

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
data.e = [];
data.xe = [];
data.length = [];

% controller parameters
ctrl.k = 2;

% show the real-time simulation

for t = 0:sim.step:sim.t
    
    if strcmp(option, 'switch')
        if t<=70
            sim.L = sim.L1;
        elseif t>140
            sim.L = sim.L3;
        else
            sim.L = sim.L2;
        end
    else
        sim.L = sim.L1;
    end
    
    sim.L = kron(sim.L, eye(2));
    
    data.cd = [];
    
    for s = 0:param.step:1
        
        if strcmp(option, 'fault')
            x = 0.01*(4*t+(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800)*cos(2*pi*s));
            y = 0.01*(2*t+(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800)*sin(2*pi*s));
        else
            x = 0.01*(4*t+(t*sin(4*pi*s)+800)*cos(2*pi*s));
            y = 0.01*(4*t+(t*cos(4*pi*s)+800)*sin(2*pi*s));
        end
        
        data.cd = [data.cd; [x; y]];
        
    end
    
    [len, len_list] = cLength(0, 1, param.step, t, option);
    param.Gs = generate_Gs_new(len_list, param.Ns);
    
    % curve fitting
    param.coff = zeros(4*param.Ns+2,1)+1;
    % param.coff = pinv(param.Gs'*param.Gs)*param.Gs'*data.cd;
    ndelta = 1e+4;
    
    index = 1;
    fs = param.Gs*param.coff;
    e = fs-data.cd;
    while ndelta>1e-4
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
    [ctrl.u, sim.coff, sim.dcoff] = controller(t, sim.N, ctrl.k, sim.L, sim.Gs, X0, ...
        param.coff, dys);
    
    % integration
    tt = [t, t+sim.step];
    
    if strcmp(dys, 'nonholonomic')
        [~, Temp] = ode45(@(t, X)dynamics2(t, X, ctrl.u, sim.N), tt, X0, ...
            odeset('RelTol', 1e-6, 'AbsTol', 1e-6));
        
        % update
        X0 = Temp(end,:)';
        sim.c = X0(1:2*sim.N);
        sim.theta = X0(2*sim.N+1:end);
        for i=1:sim.N
            sim.c(2*i-1) = sim.c(2*i-1)+0.01*cos(sim.theta(i));
            sim.c(2*i) = sim.c(2*i)+0.01*sin(sim.theta(i));
        end
        
        data.theta = [data.theta, sim.theta];
    else
        [~, Temp] = ode45(@(t, X)dynamics1(t, X, ctrl.u), tt, X0, ...
            odeset('RelTol', 1e-6, 'AbsTol', 1e-6));
        
        % update
        X0 = Temp(end,:)';
        sim.c = X0;
    end
    
    % saving data
    data.t = [data.t, t];
    data.c = [data.c, sim.c];
    data.u = [data.u, ctrl.u];
    data.dcoff = [data.dcoff, sim.dcoff];
    data.cs = param.Gs*sim.coff;
    
    data.f1.cd = [data.f1.cd, data.cd];
    data.f1.cds = [data.f1.cds, data.cds];
    data.f1.c = [data.f1.c, data.c(:,end)];
    data.f1.cs = [data.f1.cs, data.cs];
    
    cur_len = len_evol(param.G_cl, sim.coff, sim.N);
    ref_len = len_evol(param.G_cl, param.coff, sim.N);
    data.length = [data.length, [len/sim.N; cur_len; ref_len]];
    
    % plot real-time data
    %     plot(data.cd(1:2:end), data.cd(2:2:end), 'k',...
    %         data.cds(1:2:end), data.cds(2:2:end), 'r',...
    %         data.cs(1:2:end), data.cs(2:2:end),'g',...
    %         data.c(1:2:end, end), data.c(2:2:end, end), 'b*', 'LineWidth', 1);
    %     grid on;
    %     axis([-10, 25,-12, 23]);
    %     legend('Original Curve', 'Re-paramerterized Curve', 'Real-time Curve',...
    %         'Agents Real-time Positions',  'Location', 'NorthWest');
    %     drawnow;
    
    t
    
end

toc;

%% plotting
% % gif
% plotGif (data, param.Ns, option);
% % figures
% plotFigure(data, sim, option);

end
