function [u, s, ds] = controller(k, L, G, N, c, sd)
 
    s = (G'*G)\G'*c;
%     s = G'*pinv(G*G')*c;
%     s = pinv(G)*c;
    
    ds = s - sd;
    
    u = -k*L*G*ds;
%     u = -k*G*ds;
    
end

