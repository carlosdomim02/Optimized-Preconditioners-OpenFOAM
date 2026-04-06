clc; clear all;
% Reading matrix before
rows = dlmread('./MatricesResults/resultBefore/rows.txt');
columns = dlmread('./MatricesResults/resultBefore/columns.txt');
values = dlmread('./MatricesResults/resultBefore/values.txt');

matrix = sparse(rows, columns, values);

disp('Condition Number Before: ');
disp(condest(matrix));

clear all;
% Reading matrix after
rows = dlmread('./MatricesResults/resultAfterDiag/rows.txt');
columns = dlmread('./MatricesResults/resultAfterDiag/columns.txt');
values = dlmread('./MatricesResults/resultAfterDiag/values.txt');

matrix = sparse(rows, columns, values);

disp('Condition Number After Diagonal: ');
disp(condest(matrix));

clear all;
% Reading matrix after
rows = dlmread('./MatricesResults/resultAfterDIC/rows.txt');
columns = dlmread('./MatricesResults/resultAfterDIC/columns.txt');
values = dlmread('./MatricesResults/resultAfterDIC/values.txt');

matrix = sparse(rows, columns, values);

disp('Condition Number After DIC: ');
disp(condest(matrix));

%cond(matrix,1)
% si usas cond (la precisa) con una matriz sparse hace condest por defecto

% Forma precisa de hacer el mismo calculo que cond
%inv_mat = inv(matrix) % Se queda sin memoria al igual que python
%disp(norm(matrix, 1) * norm(inv_mat, 1))