function[sparse_solver_solution, end_time, A_sparse, b] = sparse_solver(nx ,ny)
start_time = tic;
A =  mat_gen(nx , ny);
b = zeros(nx ,ny);
for j = 1 : nx
    for k = 1 : ny
        b(j,k) = -2*(pi^2)*sin(pi*(j/(nx+1)))*sin(pi*(k/(ny+1)));
    end
end
b = reshape(b,nx*ny,1);
A_sparse = sparse(A);
x = A_sparse\b;
sparse_solver_solution = reshape(x , [nx ny]);
sparse_solver_solution = padarray(sparse_solver_solution , [1 1] , 0);
end_time =toc(start_time); 
end