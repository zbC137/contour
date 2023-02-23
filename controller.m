function [u, s, ds] = controller(t, N, k, L, G, X, sd, dys)

c = X(1:2*N);

%     s = (G'*G)\G'*c;
%     s = G'*pinv(G*G')*c;
s = pinv(G)*c;
ds = s - sd;
dc = c - G*sd;
uTemp = -k*L*G*ds -dc;

if strcmp(dys, 'nonholonomic')
    xl = zeros(2*N,1);
    for i=1:N
        xl(2*i-1) = X(2*i-1)+0.01*cos(X(2*N+i));
        xl(2*i) = X(2*i)+0.01*sin(X(2*N+i));
    end
    
    s = pinv(G)*xl;
    ds = s-sd;
    dc = xl-G*sd;
    u = -k*L*G*ds-dc;
    
%     theta = X(2*N+1:end);
%     u1 = zeros(N,1);
%     for i=1:N
%         u1(i) = uTemp(2*i-1)*cos(theta(i))+uTemp(2*i)*sin(theta(i));
%         if u1(i)>0.1
%             u1(i) = 0.1;
%         elseif u1(i)<-0.1
%             u1(i)=-0.1;
%         end
%     end
    
    %     u1 = u1 - dc;
%     u2 = zeros(N,1)+cos(t);
%     u = [u1; u2];
    %     u = -k*L*dc - dc;
    %     u = -k*G*ds;
else
    u = uTemp;
end

end
