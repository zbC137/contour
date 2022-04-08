function G = generate_G(step, boundary, N, flag)

G = [];

if flag==1
    for i=0:step:boundary
    
        temp = zeros(2, 4*N+2);
        for j = 1:N
        
            temp(:,4*j-3:4*j) = [cos(j*i), sin(j*i),        0,       0
                                        0,        0, cos(j*i), sin(j*i)];     
    
        end
    
        temp(:, 4*N+1:4*N+2) = eye(2);
        G = [G; temp];
    
    end
end

if flag==2
    for i=1:step
        
        temp = zeros(2, 4*N+2);
        rho = boundary(i);
        
        for j = 1:N
            temp(:,4*j-3:4*j) = [cos(j*rho), sin(j*rho),          0,          0
                                          0,          0, cos(j*rho), sin(j*rho)];
        end
        
        temp(:, 4*N+1:4*N+2) = eye(2);
        G = [G; temp];
        
    end
end

end

