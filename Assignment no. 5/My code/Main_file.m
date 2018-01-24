clc;
clear all;
close all;
step_size = [1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048, 1/4096];
N_X = [3, 7, 15, 31];
N_Y = [3, 7, 15, 31];
end_time = [1/8, 2/8, 3/8, 4/8];

%%%%%% Approximation of Instationary Heat Equation by Explicit Euler %%%%%%

for i = 1:length(end_time)
   fig = figure('name', strcat('Solution of instationary 2D-heat eqauation by Explicit Euler at time: ', num2str(end_time(i))));
    for j = 1:length(N_X)
        hx = 1/(N_X(j)+1);
        hy = 1/(N_Y(j)+1);
        for k = 1:length(step_size)
            % calculation of number of iteration to reach the perticular time
            iterations = end_time(i)/step_size(k); 
            % Making a matrix T by using initial and boundary conditions
            T = ones(N_X(j), N_Y(j));
            T = padarray(T, [1 1], 0);
            for l = 1:iterations
                T = exp_euler_2D(T, step_size(k), N_X(j), N_Y(j));
            end
            % subplot_tight is not a built-in function, we got it from here
            % https://de.mathworks.com/matlabcentral/fileexchange/30884-controllable-tight-subplot?focused=6614736&tab=function
            % It adjusts the subplots with respect to neighbouring plots
            subplot_tight(length(N_X), length(step_size), (7*(j-1))+k);
            surf(T);
            %For maximizing the plot window
            set(gcf, 'Position', get(0,'Screensize'));
            %As the plate under consideration is square with dimensions of
            %[0 1]^2 so that's why both axis are being set to 0 - 1.  
			[x_axis,y_axis] = meshgrid(0:hx:1,0:hy:1);
			surf(x_axis,y_axis,T);
            if j==1
                title(strcat('dt = 1/', num2str(1/step_size(k))));
            end
            if k==1
                zlabel(strcat('Nx=Ny = ',num2str(N_X(j))));
            end  
        end  % end time step loop
    end % end grid size loop
    saveas(fig, strcat('Plot for Explicit Euler at time = ', num2str(end_time(i)), '.jpg'));
end  % end time loop

%%%%%% Approximation of Instationary Heat Equation by Implicit Euler %%%%%%

fig2 = figure('name', 'Gauss Seidel approsimation to Heat Equation by Implicit Euler at step size 1/64');
for i = 1:length(N_X)
    dt = 1/64;
    for j = 1:length(end_time)
        % Making a matrix T by using initial and boundary conditions
        T = ones(N_X(i), N_Y(i));
        T = padarray(T, [1 1], 0);
        hx = 1/(N_X(i)+1);
        hy = 1/(N_Y(i)+1);
        % Calculations for number of steps to reach end_time(j)
        iterations = end_time(j)/dt;
        for k = 1:iterations
            T = implicit_Euler_2D(T, dt, N_X(i), N_Y(i));
        end
        subplot_tight(length(N_X), length(end_time), (4*(i-1))+j);
        surf(T);
        %For maximizing the plot window
        set(gcf, 'Position', get(0,'Screensize')); 
        %As the plate under consideration is square with dimensions of
        %[0 1]^2 so that's why both axis are being set to 0 - 1.  
        [x_axis, y_axis] = meshgrid(0:hx:1,0:hy:1);
        surf(x_axis, y_axis, T);
        if i==1
           title(strcat('time = ', num2str(j), '/8'));
        end
        if j==1
           zlabel(strcat('Nx=Ny = ',num2str(N_X(i))));
        end
    end % end of grid size loop
end % end of max time loop
saveas(fig2, strcat('Implicit Eluer Plot with dt = ', num2str(dt), '.jpg'));
% dt <= h2/4
stability_array = zeros(4,7);
for i = 1:7             %Computing the table for stability criteria
    for j = 1:4
        dt = step_size(1,i);
        h = 1/(N_X(1,j)+1);
        if (dt <= (h^2)/4)               %condition computed by Von-neumann Stability analysis
            stability_array(j,i) = 1;
        end
    end
end

T = array2table(stability_array,'VariableNames',{'dt64','dt128','dt256','dt512','dt1024','dt2048','dt4096'},'RowNames',{'N3','N7','N15','N32'})
