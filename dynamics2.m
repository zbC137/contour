function dX = dynamics2(t, X, u, N)
dX = zeros(3*N, 1);
for i=1:N
    dX(2*i-1) = u(2*i-1)*cos(X(2*N+i));
    dX(2*i) =  u(2*i)*sin(X(2*N+i));
end

dX(2*N+1:end) = u(2*N+1:end);

end
