function cur_len = len_evol(G, coff, N)

pos = G*coff;
x_pos = pos(1:2:end);
y_pos = pos(2:2:end);

len_list = sqrt((x_pos(2:end)-x_pos(1:2400)).^2+(y_pos(2:end)-y_pos(1:2400)).^2);

step = 2400/N;
cur_len = zeros(N,1);

for i=1:N
    
    cur_len(i) = sum(len_list((i-1)*step+1:i*step));
    
end

end

