F = 2*kron(sim.L2, eye(2))*sim.Gs*pinv(sim.Gs)+eye(size(sim.L, 1));
epsilon = min(eig(F));

xd = sim.Gs*data.coff;

xd1 = xd(:, 2:end);
xd2 = xd(:, 1:2000);

dxd = xd1-xd2;
a = max(vecnorm(dxd));
b = norm(data.xe(:, 1));

thre = -1/epsilon*log(1-a/b)
