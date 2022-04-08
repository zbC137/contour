function G = generate_Gs(start, step, boundary, N)

G = [];

for i=start:step:boundary
    
    temp = zeros(2, 4*N+2);
    for j = 1:N
        
        temp(:,4*j-3:4*j) =  [cos(2*pi*j*i), sin(2*pi*j*i),        0,       0
            0,        0,         cos(2*pi*j*i), sin(2*pi*j*i)];
        
    end
    
    temp(:, 4*N+1:4*N+2) = eye(2);
    G = [G; temp];
    
end
end
