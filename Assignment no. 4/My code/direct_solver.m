function [result , end_time , A , b] = direct_solver(nx , ny)
start_time = tic;
A =  mat_gen(nx , ny);
b = zeros(nx ,ny);
for j = 1 : nx
    for k = 1 : ny
        b(j,k) = -2*(pi^2)*sin(pi*(j/(nx+1)))*sin(pi*(k/(ny+1)));
    end
end
b = reshape(b,nx*ny,1);
x = A\b;
result = reshape(x , [nx ny]);
result = padarray(result , [1 1] , 0);
end_time = toc(start_time);
end