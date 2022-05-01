function [u, s, ds] = controller(k, L, G, c, sd)
 
%     s = (G'*G)\G'*c;
%     s = G'*pinv(G*G')*c;
    s = pinv(G)*c;
    
    ds = s - sd;
    
    dc = c - G*sd;
    
    u = -k*L*G*ds - dc;
%     u = -k*L*dc - dc;
%     u = -k*G*ds;
    
end

