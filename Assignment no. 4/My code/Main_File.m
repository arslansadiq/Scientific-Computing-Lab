clc;
clear all;
close all;
N_X = [7 15 31 63];
N_Y = [7 15 31 63];
computational_record_direct = zeros(2,4);
computational_record_sparse = zeros(2,4);
computational_record_gs = zeros(2,4);
for i = 1 : size(N_X , 2)
    n_x = N_X(i);
    n_y = N_Y(i);
    
%%%%%%%%%%%%%%%  Direct Solver Solution to Heat Equation   %%%%%%%%%%%%%%%%

    [direct_solution , direct_solver_time , A , b] = direct_solver(n_x , n_y);
    computational_record_direct(1 , i) = direct_solver_time;
    storage_dir_sol = whos('direct_solution');
    storage_A = whos('A');
    storage_b = whos('b');
    computational_record_direct(2 , i) = storage_dir_sol.bytes + storage_A.bytes + storage_b.bytes;
    figure(i)
    subplot(2, 1, 1);
    surf(direct_solution);
    title(strcat('3-d Plot of Direct Solver for Grid Size: ' , num2str(N_X(i))))
    hold on
    subplot(2, 1, 2);
    contour(direct_solution);
    title(strcat('Contour Plot of Direct Solver for Grid Size: ' , num2str(N_X(i))))
    hold on
    clearvars -except computational_record_direct computational_record_sparse i n_x n_y N_X N_Y computational_record_gs
    
%%%%%%%%%%%%%%%  Sparse Sover approach to Heat Equation   %%%%%%%%%%%%%%%%

    [sparse_solver_solution , Sparse_solver_time , A_Sparse , b] = sparse_solver(n_x , n_y);
    computational_record_sparse(1 , i) = Sparse_solver_time;
    storage_sparse_sol = whos('sparse_solver_solution');
    storage_A_Sparse = whos('A_Sparse');
    storage_b = whos('b');
    computational_record_sparse(2 , i) = storage_sparse_sol.bytes + storage_A_Sparse.bytes + storage_b.bytes;
    figure(i+4)
    subplot(2, 1, 1);
    surf(sparse_solver_solution);
    title(strcat('3-d Plot of Sparse Solver for Grid Size: ' , num2str(N_X(i))))
    hold on
    subplot(2, 1, 2);
    contour(sparse_solver_solution);
    title(strcat('Contour Plot of Sparse Solver for Grid Size: ' , num2str(N_X(i))))
    hold on
    clearvars -except computational_record_direct computational_record_sparse i n_x n_y N_X N_Y computational_record_gs
    
%%%%%%%%%%%%  Gauss-Seidel Solver approach to Heat Equation   %%%%%%%%%%%%%

    [gauss_seidel_sol , gs_time , b] = gauss_seidel_solver(N_X(i),N_Y(i));
    computational_record_gs(1 , i) = gs_time;
    storage_gs_sol = whos('gauss_seidel_sol');
    storage_b = whos('b');
    computational_record_gs(2 , i) = storage_gs_sol.bytes + storage_b.bytes;
    figure(i+8)
    subplot(2,1,1);
    surf(gauss_seidel_sol);
    title(strcat('3-d Plot of Gauss-Seidel Solver for Grid Size: ' , num2str(N_X(i))));
    hold on
    subplot(2, 1, 2);
    contour(gauss_seidel_sol);
    title(strcat('Contour Plot of Gauss-Seidel Solver for Grid Size: ' , num2str(N_X(i))))
    hold on
    clearvars -except computational_record_direct computational_record_sparse i n_x n_y N_X N_Y computational_record_gs
end

   %%%%%%%%%%%%%%%%%%%%%   Computational Records   %%%%%%%%%%%%%%%%%%%%

Table_Direct_Solver = array2table(computational_record_direct);
Table_Direct_Solver.Properties.RowNames = {'Run Time' , 'Storage'};
Table_Direct_Solver.Properties.VariableNames = {'Grid_Size_7' , 'Grid_Size_15' , 'Grid_Size_31' , 'Grid_Size_63'};
Table_Direct_Solver
Table_Sparse_Solver = array2table(computational_record_sparse);
Table_Sparse_Solver.Properties.RowNames = {'Run Time' , 'Storage'};
Table_Sparse_Solver.Properties.VariableNames = {'Grid_Size_7' , 'Grid_Size_15' , 'Grid_Size_31' , 'Grid_Size_63'};
Table_Sparse_Solver
Table_Gauss_Seidel_Solver = array2table(computational_record_gs);
Table_Gauss_Seidel_Solver.Properties.RowNames = {'Run Time' , 'Storage'};
Table_Gauss_Seidel_Solver.Properties.VariableNames = {'Grid_Size_7' , 'Grid_Size_15' , 'Grid_Size_31' , 'Grid_Size_63'};
Table_Gauss_Seidel_Solver

%%%%%%%%%%%%%%%%%%%%%%%    Gauss Seidel Errors    %%%%%%%%%%%%%%%%%%%%%%

N_X = [7 15 31 63 127];
N_Y = [7 15 31 63 127];
error_table = zeros(2,5);
for i=1:size(N_X,2)
    [T , ~ , ~] = gauss_seidel_solver(N_X(i) , N_Y(i));
    T_exact = zeros(N_X(i) , N_Y(i));
    hx = 1/(N_X(i)+1);
    hy = 1/(N_Y(i)+1);
    for j=1:N_X(i)
        for k=1:N_Y(i)
            T_exact(j,k) = sin(pi*j*hx)*sin(pi*k*hy);
        end
    end
    error = sqrt(1/(N_X(i)*N_Y(i)) * sum(sum((T_exact - T(2:N_X(i)+1 , 2:N_Y(i)+1)).^2)));
    error_table(1,i) = error;
    if i ~= 1
        error_table(2,i) = error_table(1,i-1)/error_table(1,i);
    end
end
Gauss_Seidel_Error_Table = array2table(error_table);
Gauss_Seidel_Error_Table.Properties.RowNames = {'Error' , 'Error Reduced'};
Gauss_Seidel_Error_Table.Properties.VariableNames = {'Grid_Size_7', 'Grid_Size_15', 'Grid_Size_31', 'Grid_Size_63',...
    'Grid_Size_127'};
Gauss_Seidel_Error_Table