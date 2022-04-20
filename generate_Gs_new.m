function G = generate_Gs_new(len_list, N)

G = [];
num = size(len_list, 2);

for i=1:num-1
    
    l = sum(len_list(1:i));
    
    Gs_temp = zeros(2, 4*N+2);
    for j=1:N
        
        Gs_temp(:,4*j-3:4*j) = [cos(2*pi*j*l), sin(2*pi*j*l),        0,       0
                                           0,        0,         cos(2*pi*j*l), sin(2*pi*j*l)];
                                    
    end
    
    Gs_temp(:, 4*N+1:4*N+2) = eye(2);
    G = [G; Gs_temp];
    
end
end
