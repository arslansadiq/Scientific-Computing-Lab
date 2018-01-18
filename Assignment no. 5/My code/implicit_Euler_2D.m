function Tout = implicit_Euler_2D(X, dt, nx, ny)
hx = 1/(nx+1);
hy = 1/(ny+1);
Tout = zeros(nx+2, ny+2);
a = dt/hx^2;
b = dt/hy^2;
c = 1 + 2*(dt/hx^2) + 2*(dt/hy^2);
residual = 5;
N = nx*ny;
while residual>10^-6
    %Iteration for Gauss Seidel 
    % Equation for Implicit Euler:
    % Tout = X + dt*Tout
    for i = 2:nx+1
        for j = 2:ny+1
            Tout(i, j) = (a/c)*(Tout(i+1, j) + Tout(i-1, j))+...
                    (b/c)*(Tout(i, j+1) + Tout(i, j-1))+...
                    (1/c)*(X(i,j));
        end
    end
    % Calculation residual at every point. It Should as close to zero as
    % possible
    error = zeros(nx, ny);
    for i = 2:nx+1
        for j=2:ny+1
            error(i-1, j-1) = X(i, j) + (a)*(Tout(i+1, j) + Tout(i-1, j))+...
                (b)*(Tout(i, j+1) + Tout(i, j-1))-...
                (c)*(Tout(i, j)); 
        end
    end
    residual = sqrt((1/N)*(sum(sum(error.^2))));
end
end