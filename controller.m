function [u, s, ds] = controller(t, N, k, L, G, X, sd, dys)

c = X(1:2*N);

%     s = (G'*G)\G'*c;
%     s = G'*pinv(G*G')*c;
s = pinv(G)*c;
ds = s - sd;
dc = c - G*sd;

if strcmp(dys, 'nonholonomic')
    theta = X(2*N+1:end);
    u1 = -k*L*G*ds -dc;
    for i=1:N
        u1(2*i-1) = u1(2*i-1)*cos(theta(i));
        u1(2*i) = u1(2*i)*sin(theta(i));
    end
    
    %     u1 = u1 - dc;
    u2 = zeros(N,1)+cos(t);
    u = [u1; u2];
    %     u = -k*L*dc - dc;
    %     u = -k*G*ds;
else
    u = -k*L*G*ds -dc;
end

end
