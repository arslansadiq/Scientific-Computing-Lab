function A = mat_gen(N_X , N_Y)
A = zeros(N_X*N_Y);
h_x = 1/(1 + N_X);
h_y = 1/(1 + N_Y);
for j = 1 : N_X*N_Y
    A(j,j) = -2*((1/h_x)^2 + (1/h_y^2));
end
for j = 1 : (N_X-1)
    for k = 1 : (N_Y)
        A(j+(k-1)*N_X , j+(k-1)*N_X+1) = 1*(1/h_x)^2;
        A(j+(k-1)*N_X+1 , j+(k-1)*N_X) = 1*(1/h_x)^2;
    end
end
for i = 1 : N_X
    for j = 1 : (N_Y-1)
        A(i+(j-1)*N_X , i+j*N_X) = 1*(1/h_y)^2;
        A(i+j*N_X , i+(j-1)*N_X) = 1*(1/h_y)^2;
    end
end
end