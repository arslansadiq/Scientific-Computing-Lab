clc;
close all;
clear all;
time_end = 5;
time = 0 : 0.001 : time_end;
time_step = [1 , 0.5 , 0.25 , 0.125];
exact_sol = 10./(1+9*exp(-time));

%%%%%%%%%%%          Implementation of Euler Method       %%%%%%%%%%%%%
y0 = 1;
exact_error_euler = zeros(1 , size(time_step , 2));
approx_error_euler = zeros(1 , size(time_step , 2));
best_euler_approx = euler(time_step(4) , time_end , y0);
for i=1:size(time_step , 2)
    i
    tic
    euler_approx = euler(time_step(i) ,time_end , y0);
    toc
    figure(2)
    plot(0 : time_step(i) : time_end , euler_approx , 'DisplayName' , strcat('Time Step= ' , string(time_step(i))));
    title('Numerical Method: Euler Method Approximation')
    hold on;
    exact_solution_E = 10./(1+9*exp(-(0 : time_step(i) : time_end)));
    exact_error_euler(i) = error(euler_approx , time_step(i) , time_end , exact_solution_E);
    approx_error_euler(i) = approx_error(euler_approx , time_step(i) , time_end , best_euler_approx , i);
end
plot(time , exact_sol , 'DisplayName' , 'Analytical Solution');
legend show
legend('Location','northwest')
error_reduced_euler = zeros(1 , size(time_step , 2));
error_reduced_euler(1) = 0.0;
for i = 2:size(error_reduced_euler , 2)
    error_reduced_euler(i) = exact_error_euler(i-1)/exact_error_euler(i);
end
table_vector_euler = [time_step; exact_error_euler; error_reduced_euler; approx_error_euler];
Euler_Table = array2table(table_vector_euler , 'RowNames' , {'dt','error','error red','error app'})

%%%%%%%%%%%          Implementation of Heun Method       %%%%%%%%%%%%%
exact_error_heun = zeros(1 , size(time_step , 2));
approx_error_heun = zeros(1 , size(time_step , 2));
best_heun_approx = Heun(time_step(4) , time_end , y0);
for i=1:size(time_step , 2)
    i
    tic
    heun_approx = Heun(time_step(i) , time_end , y0);
    toc
    figure(3)
    plot(0 : time_step(i) : time_end , heun_approx , 'DisplayName' , strcat('Time Step= ' , string(time_step(i))));
    title('Numerical Method: Heun Method Approximation')
    hold on;
    exact_solution_H = 10./(1+9*exp(-(0 : time_step(i) : time_end)));
    exact_error_heun(i) = error(heun_approx , time_step(i) , time_end , exact_solution_H);
    approx_error_heun(i) = approx_error(heun_approx , time_step(i) , time_end , best_heun_approx , i);
end
plot(time , exact_sol , 'DisplayName' , 'Analytical Solution');
legend show
legend('Location','northwest')
error_reduced_heun = zeros(1 , size(time_step , 2));
error_reduced_heun(1) = 0.0;
for i = 2:size(error_reduced_heun , 2)
    error_reduced_heun(i) = exact_error_heun(i-1)/exact_error_heun(i);
end
table_vector_heun = [time_step; exact_error_heun; error_reduced_heun; approx_error_heun];
Heun_Table = array2table(table_vector_heun , 'RowNames' , {'dt','error','error red','error app'})

%%%%%%%%%%%          Implementation of Runge-Kutta Method       %%%%%%%%%%%%%
exact_error_rk = zeros(1 , size(time_step , 2));
approx_error_rk = zeros(1 , size(time_step , 2));
best_runge_kutta_approx = Runge_Kutta(time_step(4) , time_end , y0);
for i=1:size(time_step , 2)
    i
    tic
    runge_kutta_approx = Runge_Kutta(time_step(i) , time_end , y0);
    toc
    figure(4)
    plot(0 : time_step(i) : time_end , runge_kutta_approx , 'DisplayName' , strcat('Time Step= ' , string(time_step(i))));
    title('Numerical Method: Runge-Kutta Method Approximation')
    hold on;
    exact_solution_rk = 10./(1+9*exp(-(0 : time_step(i) : time_end)));
    exact_error_rk(i) = error(runge_kutta_approx , time_step(i) , time_end , exact_solution_rk);
    approx_error_rk(i) = approx_error(runge_kutta_approx , time_step(i) , time_end , best_runge_kutta_approx , i);
end
plot(time , exact_sol , 'DisplayName' , 'Analytical Solution');
legend show
legend('Location','northwest')
error_reduced_rk = zeros(1 , size(time_step , 2));
error_reduced_rk(1) = 0.0;
for i = 2:size(error_reduced_rk , 2)
    error_reduced_rk(i) = exact_error_rk(i-1)/exact_error_rk(i);
end
table_vector_rk = [time_step; exact_error_rk; error_reduced_rk; approx_error_rk];
RK_Table = array2table(table_vector_rk , 'RowNames' , {'dt','error','error red','error app'})

   %%%%%%%%%%%          Euler Method's Definition        %%%%%%%%%%%%%

function approx_value = euler(dt , time_end , y0)
approx_value =  zeros(1 , (time_end./dt)+1);
approx_value(1) = y0;
for i = 1:(size(approx_value , 2)-1)
    approx_value(i+1) = approx_value(i) + (dt)*((1 - (approx_value(i)/10))*approx_value(i));
end
end

   %%%%%%%%%%%          Heun Method's Definition        %%%%%%%%%%%%%

function approx_value = Heun(dt , time_end , y0)
approx_value = zeros(1 , (time_end./dt)+1);
approx_value(1) = y0;
for i = 1:(size(approx_value , 2)-1)
    y_derivative_at_i = ((1 - ((approx_value(i)/10)))*approx_value(i));
    y_tmp = approx_value(i) + (dt)*(y_derivative_at_i);
    y_derivative_at_i_plus_1 = ((1 - (y_tmp/10))*y_tmp);
    approx_value(i+1) = approx_value(i) + (dt)*(0.5)*(y_derivative_at_i + y_derivative_at_i_plus_1);
end
end

   %%%%%%%%%%%          Runge Kutta Method's Definition        %%%%%%%%%%%%%
function approx_value = Runge_Kutta(dt , time_end , y0)
approx_value = zeros(1 , (time_end./dt)+1);
approx_value(1) = y0;
for i = 1:(size(approx_value , 2)-1)
    Y1 = (1 - (approx_value(i)/10))*approx_value(i);        %Y1
    y2_tmp = approx_value(i) + ((dt/2)*(Y1));
    Y2 = (1 - (y2_tmp/10))*y2_tmp;                          %Y2
    y3_tmp = approx_value(i) + ((dt/2)*(Y2));
    Y3 = (1 - (y3_tmp/10))*y3_tmp;                          %Y3
    y4_tmp = approx_value(i) + ((dt)*(Y3));
    Y4 = (1 - (y4_tmp/10))*y4_tmp;
    approx_value(i+1) = approx_value(i) + ((dt)*(1/6)*(Y1 + 2*Y2 + 2*Y3 + Y4));
end
end

%%%%%%%%%%          Function for Exact Error     %%%%%%%%%%%%%%%

function e_error = error(approximation_vector , dt , time_end , exact_solution)
y = sum((approximation_vector - exact_solution).^2);
e_error = sqrt((dt/time_end)*y);
end

%%%%%%%%%%        Function for the computation of approximation error     %%%%%%%%%%%%%%%

function a_error = approx_error(approximation_vector , dt , time_end , best_solution , i)
y = sum((approximation_vector - downsample(best_solution , 2.^(4-i))).^2);
a_error = sqrt((dt/time_end)*y);
end
