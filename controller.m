function [u, s, ds] = controller(k, L, G, N, c, sd)

%     if size(G, 1)>N
%         s = (G'*G)\G'*c;
%     elseif size(G, 1)<N
%         s = G'*((G*G')\c);
%     else
%         s = G\c;
%     end    
%     s = pinv(G'*G)*G'*c;
    s = G'*pinv(G*G')*c;
% s = pinv(G)*c;
    
    ds = s - sd;
    
    u = -k*L*G*ds;
%     u = -k*G*ds;
    
end

