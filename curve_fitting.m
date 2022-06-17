clc; clear;

%% parameters
param.H = [6, 12, 18];
param.N = 15;

param.t = 200;
param.dt = 0.1;

param.step = 0.001;

%% curve fitting
data.t = 0:param.dt:param.t;
data.length1 = [];
data.length2 = [];

figure(1)

for i = 1:3
    i
    
    h = param.H(i);
    param.G_cl = generate_Gs(0, 1/2400, 1, h);
    
    for t = 0:param.dt:param.t
        
        t
        
        data.cd1 = [];
        data.cd2 = [];
        
        for s = 0:param.step:1
            
            x1 = 0.01*(4*t+(t*sin(4*pi*s)+800)*cos(2*pi*s));
            y1 = 0.01*(4*t+(t*cos(4*pi*s)+800)*sin(2*pi*s));
            data.cd1 = [data.cd1; [x1; y1]];
            
            x2 = 0.01*(4*t+(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800)*cos(2*pi*s));
            y2 = 0.01*(2*t+(t*sin(4*pi*s)+2*t*cos(10*pi*s)+800)*sin(2*pi*s));
            data.cd2 = [data.cd2; [x2; y2]];
            
        end
        
        [len1, len_list1] = cLength(0, 1, param.step, t, 'sim');
        param.Gs1 = generate_Gs_new(len_list1, h);
        
        [len2, len_list2] = cLength(0, 1, param.step, t, 'fault');
        param.Gs2 = generate_Gs_new(len_list2, h);
        
        % curve fitting
        param.coff1 = zeros(4*h+2,1)+1;
        param.coff2 = zeros(4*h+2,1)+1;
        ndelta1 = 1e+4;
        ndelta2 = 1e+4;
        
        index = 1;
        fs1 = param.Gs1*param.coff1;
        e1 = fs1-data.cd1;
        while ndelta1>1e-4
            
            index;
            
            delta = -pinv(param.Gs1'*param.Gs1)*param.Gs1'*e1;
            param.coff1 = param.coff1+delta;
            fs1 = param.Gs1*param.coff1;
            e1 = fs1-data.cd1;
            error = norm(e1);
            ndelta1 = norm(delta);
            
            index = index+1;
            
        end
        
        data.cds1 = param.Gs1*param.coff1;
        
        index = 1;
        fs2 = param.Gs2*param.coff2;
        e2 = fs2-data.cd2;
        while ndelta2>1e-4
            
            index;
            
            delta = -pinv(param.Gs2'*param.Gs2)*param.Gs2'*e2;
            param.coff2 = param.coff2+delta;
            fs2 = param.Gs2*param.coff2;
            e2 = fs2-data.cd2;
            error = norm(e2);
            ndelta2 = norm(delta);
            
            index = index+1;
            
        end
        data.cds2 = param.Gs2*param.coff2;
        
        ref_len1 = len_evol(param.G_cl, param.coff1, param.N);
        data.length1 = [data.length1, [len1/param.N; ref_len1]];
        
        ref_len2 = len_evol(param.G_cl, param.coff2, param.N);
        data.length2 = [data.length2, [len2/param.N; ref_len2]];
        
        % plot real-time data   
        subplot(1,2,1);
        plot(data.cd1(1:2:end), data.cd1(2:2:end), 'k',...
            data.cds1(1:2:end), data.cds1(2:2:end), 'r', 'LineWidth', 1);
        grid on;
        axis([-10, 20, -10, 15]);
        
        subplot(1,2,2);
        plot(data.cd2(1:2:end), data.cd2(2:2:end), 'k',...
            data.cds2(1:2:end), data.cds2(2:2:end), 'r', 'LineWidth', 1);
        grid on;
        axis([-10, 25, -10, 25]);
        
        drawnow;
            
    end
end

%% data processing and plotting
for j = 1:param.N
    
    data.length1(j+1,:) = data.length1(j+1,:) - data.length1(1,:);
    data.length2(j+1,:) = data.length2(j+1,:) - data.length2(1,:);
    
end

temp1 = data.length1(2:end, :);
temp2 = data.length2(2:end, :);

nrsme1 = sqrt(sum(temp1.^2)/param.N)./data.length1(1,:);
nrsme2 = sqrt(sum(temp2.^2)/param.N)./data.length2(1,:);

num = param.t/param.dt+1;
plt.rsme1 = [nrsme1(1,1:num); nrsme1(1, num+1:2*num); nrsme1(1, 2*num+1:end)];
plt.rsme2 = [nrsme2(1,1:num); nrsme2(1, num+1:2*num); nrsme2(1, 2*num+1:end)];
plt.t = data.t;

figure(2)
plot(plt.t, plt.rsme1(:,:), 'LineWidth', 1);
hold on;
plot(plt.t, plt.rsme2(:,:), 'LineWidth', 1);
grid on;
