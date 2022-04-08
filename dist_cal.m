function dist = dist_cal(start, final, step, N, s)

dist = 0;

for i=start:step:final
    
    G_dot = zeros(2, 4*N);
    for j = 1:N
        
       G_dot(:,4*j-3:4*j) = [-sin(j*i), cos(j*i),         0,       0
                                    0,        0, -sin(j*i), cos(j*i)];
       G_dot(:,4*j-3:4*j) = j*G_dot(:,4*j-3:4*j);
    
    end
    
    pos = G_dot*s;
    
    dist = dist + sqrt(pos(1)^2 + pos(2)^2);
    
end

dist = step*dist;

if start==final
    dist =0;
end

end

