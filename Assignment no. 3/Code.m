clc;
close all;
clear all;
time_end = 5;
time = 1 : 0.001 : time_end;
initial_condition = 20;
time_step = [0.5, 0.25, 0.125, 0.0625, 0.03125];
analytical_solution = 200./(20 - 10*exp(-7*time));
stability_table = [time_step;ones(6,5)]; 

%%%%%%%%%%%%       Implementation of Explicit Euler Method's from Worksheet 2      %%%%%%%%%%%%%

exact_error_EE = zeros(1 , size(time_step , 2));
for i=1:size(time_step , 2)
    euler_exp_approx = explicit_euler(time_step(i) , time_end , initial_condition);
    figure(1)
    axis([0 5 0 20]);
    title('Numerical Approximation Method: Explicit Euler Method');
    plot(0:time_step(i):time_end , euler_exp_approx , 'DisplayName' , strcat('Time Step: ' , string(i)));
    exact_solution_EE = 200./(20-10*exp(-7*(0 : time_step(i) : time_end)));
    exact_error_EE(i) = Error(euler_exp_approx , time_step(i) , time_end , exact_solution_EE);
    if exact_error_EE(i)>2
        stability_table(2,i) = 0;
    end
    hold on;
end
plot(time , analytical_solution , 'DisplayName' , 'Analytical Solution');
legend show
Explicit_Euler_Table = array2table(Error_Table(exact_error_EE) , 'RowNames' , {'dt','error','error red'})

%%%%%%%%%%%%       Implementation of Heun Method's from Worksheet 2      %%%%%%%%%%%%%

exact_error_EH = zeros(1 , size(time_step , 2));
for i=1:size(time_step , 2)
    heun_exp_approx = heun(time_step(i) , time_end , initial_condition);
    figure(2)
    axis([0 5 0 20]);
    title('Numerical Approximation Method: Heun Method');
    plot(0:time_step(i):time_end , heun_exp_approx , 'DisplayName' , strcat('Time Step: ' , string(i)));
    exact_solution_EH = 200./(20-10*exp(-7*(0 : time_step(i) : time_end)));
    exact_error_EH(i) = Error(heun_exp_approx , time_step(i) , time_end , exact_solution_EH);
    if exact_error_EH(i)>2
        stability_table(3,i) = 0;
    end
    hold on;
end
plot(time , analytical_solution , 'DisplayName' , 'Analytical Solution');
legend show
Heun_Table = array2table(Error_Table(exact_error_EH) , 'RowNames' , {'dt','error','error red'})

%%%%%%%%%%%%      Implementation for Implicit Euler Method's        %%%%%%%%%%%%

exact_error_IE = zeros(1 , size(time_step , 2));
for i=1:size(time_step , 2)
    euler_imp_approx = implicit_eluer(time_step(i) , time_end , initial_condition);
    figure(3)
    axis([0 5 0 20]);
    title('Numerical Approximation Method: Implicit Euler Method');
    plot(0:time_step(i):time_end , euler_imp_approx , 'DisplayName' , strcat('Time Step: ' , string(i)));
    exact_solution_IE = 200./(20-10*exp(-7*(0 : time_step(i) : time_end)));
    exact_error_IE(i) = Error(euler_imp_approx , time_step(i) , time_end , exact_solution_IE);
    if exact_error_IE(i)>2
        stability_table(4,i) = 0;
    end
    hold on;
end
plot(time , analytical_solution , 'DisplayName' , 'Analytical Solution');
legend show
Implicit_Euler_Table = array2table(Error_Table(exact_error_IE) , 'RowNames' , {'dt','error','error red'})

%%%%%%%%%%%   Implementation of 2nd Order Adam Moulton's Method's Definition   %%%%%%%%%%%%%

exact_error_adam_moltn = zeros(1 , size(time_step , 2));
for i=1:size(time_step , 2)
    adam_moulton_approx = adam_moulton(time_step(i) , time_end , initial_condition);
    figure(4)
    axis([0 5 0 20]);
    title('Numerical Approximation Method: 2nd Order Adam Moulton Method');
    plot(0:time_step(i):time_end , adam_moulton_approx , 'DisplayName' , strcat('Time Step: ' , string(i)));
    exact_solution_AM = 200./(20-10*exp(-7*(0 : time_step(i) : time_end)));
    exact_error_adam_moltn(i) = Error(adam_moulton_approx , time_step(i) , time_end , exact_solution_AM);
    if exact_error_adam_moltn(i)>2 || isnan(exact_error_adam_moltn(i))
        stability_table(5,i) = 0;
    end
    hold on;
end
plot(time , analytical_solution , 'DisplayName' , 'Analytical Solution');
legend show
Adam_Moulton_Table = array2table(Error_Table(exact_error_adam_moltn) , 'RowNames' , {'dt','error','error red'})

%%%%%%%%%%    Implementation of Linearised Version (1) of Adam Moulton's Method    %%%%%%%%%%%%%%

exact_error_adam_moltn_L1 = zeros(1 , size(time_step , 2));
for i=1:size(time_step , 2)
    L1_adam_moulton_approx = L1_adam_moulton(time_step(i) , time_end , initial_condition);
    figure(5)
    axis([0 5 0 20]);
    title('Numerical Approximation for Linearised Version (1) of Adam Moulton Method');
    plot(0:time_step(i):time_end , L1_adam_moulton_approx , 'DisplayName' , strcat('Time Step: ' , string(i)));
    exact_solution_AM_L1 = 200./(20-10*exp(-7*(0 : time_step(i) : time_end)));
    exact_error_adam_moltn_L1(i) = Error(L1_adam_moulton_approx , time_step(i) , time_end , exact_solution_AM_L1);
    if exact_error_adam_moltn_L1(i)>2
        stability_table(6,i) = 0;
    end
    hold on;
end
plot(time , analytical_solution , 'DisplayName' , 'Analytical Solution');
legend show
Adam_Moulton_L1_Table = array2table(Error_Table(exact_error_adam_moltn_L1) , 'RowNames' , {'dt','error','error red'})

%%%%%%%%%%    Implementation of Linearised Version (2) of Adam Moulton's Method    %%%%%%%%%%%%%%

exact_error_adam_moltn_L2 = zeros(1 , size(time_step , 2));
for i=1:size(time_step , 2)
    L2_adam_moulton_approx = L2_adam_moulton(time_step(i) , time_end , initial_condition);
    figure(6)
    axis([0 5 -40 40]);
    title('Numerical Approximation for Linearised Version (2) of Adam Moulton Method');
    plot(0:time_step(i):time_end , L2_adam_moulton_approx , 'DisplayName' , strcat('Time Step: ' , string(i)));
    hold on;
    exact_solution_AM_L2 = 200./(20-10*exp(-7*(0 : time_step(i) : time_end)));
    exact_error_adam_moltn_L2(i) = Error(L2_adam_moulton_approx , time_step(i) , time_end , exact_solution_AM_L2);
    if exact_error_adam_moltn_L2(i)>2
        stability_table(7,i) = 0;
    end
end
plot(time , analytical_solution , 'DisplayName' , 'Analytical Solution');
legend show
legend('Location' , 'southwest')
Adam_Moulton_L2_Table = array2table(Error_Table(exact_error_adam_moltn_L2) , 'RowNames' , {'dt','error','error red'})
Stability_Table = array2table(stability_table , 'RowNames' , {'Time Step','Explicit Euler','Heun','Implicit Euler',...
    'Adams Moulton','Adams Moulton L1','Adams Moulton L12'})
disp('Following are the pointers for Stability table');
disp('0 -> Unstability');
disp('1 -> Stability');
%%%%%%%%%%%          Explicit Euler Method's Definition        %%%%%%%%%%%%%

function approx_value = explicit_euler(dt , time_end , y0)
approx_value =  zeros(1 , (time_end./dt)+1);
approx_value(1) = y0;
for i = 1:(size(approx_value , 2)-1)
    approx_value(i+1) = approx_value(i) + (dt)*(7*(1 - (approx_value(i)/10))*approx_value(i));
end
end

%%%%%%%%%%%          Heun Method's Definition        %%%%%%%%%%%%%

function approx_value = heun(dt , time_end , y0)
approx_value = zeros(1 , (time_end./dt)+1);
approx_value(1) = y0;
for i = 1:(size(approx_value , 2)-1)
    y_derivative_at_i = (7*(1 - (approx_value(i)/10))*approx_value(i));
    y_tmp = approx_value(i) + (dt)*(y_derivative_at_i);
    y_derivative_at_i_plus_1 = (7*(1 - (y_tmp/10))*y_tmp);
    approx_value(i+1) = approx_value(i) + (dt)*(0.5)*(y_derivative_at_i + y_derivative_at_i_plus_1);
end
end

%%%%%%%%%%%          Implicit Euler Method's Definition        %%%%%%%%%%%%%

function approx_value = implicit_eluer(dt , time_end , y0)
approx_value =  zeros(1 , (time_end./dt)+1);
approx_value(1) = y0;
for i = 1:(size(approx_value , 2)-1)
    approx_value(i+1) = newton_raphson_method(10^-4 , approx_value(i) , dt , @Gx , @Gxp);
end
end

%%%%%%%%%%%          2nd Order Adam Moulton Method's Definition        %%%%%%%%%%%%%

function approx_value = adam_moulton(dt , time_end , y0, a)
approx_value =  zeros(1 , (time_end./dt)+1);
approx_value(1) = y0;
for i = 1:(size(approx_value , 2)-1)
    approx_value(i+1) = newton_raphson_method(10^-4 , approx_value(i) , dt , @Fx , @Fxp);
end
end

%%%%%%%%%%%          Linearised (1) Adam Moulton Method's Definition        %%%%%%%%%%%%%

function approx_value = L1_adam_moulton(dt , time_end , y0)
approx_value =  zeros(1 , (time_end./dt)+1);
approx_value(1) = y0;
for i = 1:(size(approx_value , 2)-1)
    d = 1 + 0.35*dt*approx_value(i);
    approx_value(i+1) = (approx_value(i) + dt*(3.5*(1 - approx_value(i)/10)*approx_value(i)) + 3.5*(dt)*approx_value(i))/d;
end
end

%%%%%%%%%%%          Linearised (1) Adam Moulton Method's Definition        %%%%%%%%%%%%%

function approx_value = L2_adam_moulton(dt , time_end , y0)
approx_value =  zeros(1 , (time_end./dt)+1);
approx_value(1) = y0;
for i = 1:(size(approx_value , 2)-1)
    d = 1 + 0.35*dt*approx_value(i) - 3.5*dt;
    approx_value(i+1) = (approx_value(i) + dt*(3.5*(1 - approx_value(i)/10)*approx_value(i)))/d;
end
end

%%%%%%%%%%%          Newton Raphson Method's Definition        %%%%%%%%%%%%%

function yi = newton_raphson_method(accuracy_limit , yi_1 , h , g , gp)
error = 5;
initial_guess = yi_1+25*rand;
itirations=0;
while (error >= accuracy_limit) && itirations < 300
    initial_guess = initial_guess - g(initial_guess , yi_1 , h)/gp(initial_guess , h);
    error = abs(g(initial_guess , yi_1 , h));
    itirations=itirations+1;
end
if error>accuracy_limit
    yi = NaN;
else
    yi = initial_guess;
end
end

%%%%%%%%%%%      Function for G(x) for Implicit Euler's Method      %%%%%%%%%%%

function g = Gx(guess , y , curr_step)
g = 0.7*curr_step*guess^2 + (1-7*curr_step)*guess - y;
end

%%%%%%%%%%%      Function for G'(x) for Implicit Euler's Method      %%%%%%%%%%%

function gp = Gxp(guess , curr_step)
gp = 1.4*guess*curr_step + 1 - 7*curr_step;
end

%%%%%%%%%%%      Function for G(x) for Adam Moulton's Method      %%%%%%%%%%%

function g = Fx(guess , y , curr_step)
g = -guess + y + ((curr_step/2)*(7*y - 0.7*y^2 + 7*guess - 0.7*guess^2));
end

%%%%%%%%%%%      Function for G'(x) for Adam Moulton's Method      %%%%%%%%%%%

function g = Fxp(guess , curr_step)
g = -1 + curr_step/2*(7 - 1.4*guess);
end

%%%%%%%%%%          Function for Exact Error     %%%%%%%%%%%%%%%

function e_error = Error(approximation_vector , dt , time_end , exact_solution)
y = sum((approximation_vector - exact_solution).^2);
e_error = sqrt((dt/time_end)*y);
end

%%%%%%%%%          Function for Creating Tables    %%%%%%%%%%%%%%%

function table = Error_Table(exact_error)
error_reduced = zeros(1 , size([0.5, 0.25, 0.125, 0.0625, 0.03125] , 2));
error_reduced(1) = 0.0;
for i = 2:size(error_reduced , 2)
    error_reduced(i) = exact_error(i-1)/exact_error(i);
end
table = [[0.5, 0.25, 0.125, 0.0625, 0.03125]; exact_error; error_reduced];
end
