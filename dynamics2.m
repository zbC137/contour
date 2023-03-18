function dX = dynamics2(t, X, u, N)

dX = zeros(3*N, 1);
% for i=1:N
%     dX(2*i-1) = u(i)*cos(X(2*N+i));
%     dX(2*i) =  u(i)*sin(X(2*N+i));
% end
% 
% dX(2*N+1:end) = u(N+1:end);
for i=1:N
    R = [cos(X(2*N+i)), -0.01*sin(X(2*N+i)); sin(X(2*N+i)), 0.01*cos(X(2*N+i))];
    ui = R\u(2*i-1:2*i);
%     if ui(1)>0.2
%         ui(1) = 0.2;
%     elseif ui(1)<-0.2
%         ui(1)=-0.2;
%     end
%     
    dX(2*i-1) = ui(1)*cos(X(2*N+i));
    dX(2*i) = ui(1)*sin(X(2*N+i));
    dX(2*N+i) = ui(2);
end

end
