function[T , end_time , b] = gauss_seidel_solver(nx , ny)
start_time = tic;
hx = 1/(1+nx);
hy = 1/(1+ny);
b = zeros(nx ,ny);
for j = 1 : nx
    for k = 1 : ny
        b(j,k) = -2*(pi^2)*sin(pi*(j/(nx+1)))*sin(pi*(k/(ny+1)));
    end
end
b_gs = padarray(b, [1 1], 0);
b = reshape(b,nx*ny,1);
T = zeros(nx+2 , ny+2);
residual = 5;
while (residual > 10^-4)
    for i = 2 : nx+1
        for j = 2 : ny+1
            T(i,j) = (((T(i-1,j) + T(i+1,j))*(1/hx)^2 + (T(i,j-1) +T(i,j+1))*(1/hy)^2) - b_gs(i,j))/(2*((1/hx)^2+(1/hy)^2));
        end
    end
    error = zeros(nx,ny);
    for i = 2 : nx+1
        for j = 2 : ny+1
            error(i-1 , j-1) = (T(i-1,j) + T(i+1,j))*(1/hx)^2 + (T(i,j-1) +T(i,j+1))*(1/hy)^2 - 2*((1/hx)^2+(1/hy)^2)*T(i,j);
        end
    end
    error = reshape(error,nx*ny,1);
    residual = sqrt(1/(nx*ny) * sum((b-error).^2));
end
end_time = toc(start_time);
end