a1 = [data.f1.cd(1:2:end,1), data.f1.cd(2:2:end,1)];
a2 = [data.f1.cds(1:2:end,1), data.f1.cds(2:2:end,1)];
a3 = [data.f1.cs(1:2:end,1), data.f1.cs(2:2:end,1)];
a4 = [data.c0(1:2:end), data.c0(2:2:end)];
a5 = data.theta(:, 1);

b1 = [data.f1.cd(1:2:end,701), data.f1.cd(2:2:end,701)];
b2 = [data.f1.cds(1:2:end,701), data.f1.cds(2:2:end,701)];
b3 = [data.f1.cs(1:2:end,701), data.f1.cs(2:2:end,701)];
b4 = [data.f1.c(1:2:end, 701), data.f1.c(2:2:end, 701)];
b5 = data.theta(:, 701);

c1= [data.f1.cd(1:2:end,1401), data.f1.cd(2:2:end,1401)];
c2 = [data.f1.cds(1:2:end,1401), data.f1.cds(2:2:end,1401)];
c3 = [data.f1.cs(1:2:end,1401), data.f1.cs(2:2:end,1401)];
c4 = [data.f1.c(1:2:end, 1401), data.f1.c(2:2:end, 1401)];
c5 = data.theta(:, 1401);

d1 = [data.f1.cd(1:2:end,end), data.f1.cd(2:2:end,end)];
d2 = [data.f1.cds(1:2:end,end), data.f1.cds(2:2:end,end)];
d3 = [ data.f1.cs(1:2:end,end), data.f1.cs(2:2:end,end)];
d4 = [data.f1.c(1:2:end, end), data.f1.c(2:2:end, end)];
d5 = data.theta(:, end);

t = data.t';
dcoff = abs(data.dcoff');
e = abs(data.e');
