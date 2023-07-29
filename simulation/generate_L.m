function l = generate_L(N, k)

if k == 0
    l = eye(N);
else
    l = eye(N)*k*2;
    
    for i = 1:N
        for j = i-k:i+k
            if j==i
                continue;
            elseif j<=0
                temp = N + j;
            elseif j>N
                temp = j - N;
            else
                temp = j;
            end
            
            l(i, temp) = -1;
        end
    end
end
end

