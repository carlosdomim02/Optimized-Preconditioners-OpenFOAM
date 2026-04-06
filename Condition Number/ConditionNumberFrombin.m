clc; clear all;
% Reading matrix before
fileRows = fopen('./MatricesResults/resultBefore/rows.bin', 'rb');
fileColumns = fopen('./MatricesResults/resultBefore/columns.bin', 'rb');
fileValues = fopen('./MatricesResults/resultBefore/values.bin', 'rb');
rows = fread(fileRows, 'int32');
columns = fread(fileColumns, 'int32');
values = fread(fileValues, 'double');
fclose(fileRows);  
fclose(fileColumns);
fclose(fileValues);

matrix = sparse(rows, columns, values);

disp('Condition Number Before: ');
disp(condest(matrix));

clear all;
% Reading matrix after
fileRows = fopen('./MatricesResults/resultAfterDIC/rows.bin', 'rb');
fileColumns = fopen('./MatricesResults/resultAfterDIC/columns.bin', 'rb');
fileValues = fopen('./MatricesResults/resultAfterDIC/values.bin', 'rb');
rows = fread(fileRows, 'int32');
columns = fread(fileColumns, 'int32');
values = fread(fileValues, 'double');
fclose(fileRows);
fclose(fileColumns);
fclose(fileValues);

matrix = sparse(rows, columns, values);

disp('Condition Number After DIC: ');
disp(condest(matrix));

clear all;
% Reading matrix after
fileRows = fopen('./MatricesResults/resultAfterDiag/rows.bin', 'rb');
fileColumns = fopen('./MatricesResults/resultAfterDiag/columns.bin', 'rb');
fileValues = fopen('./MatricesResults/resultAfterDiag/values.bin', 'rb');
rows = fread(fileRows, 'int32');
columns = fread(fileColumns, 'int32');
values = fread(fileValues, 'double');
fclose(fileRows);
fclose(fileColumns);
fclose(fileValues);

matrix = sparse(rows, columns, values);

disp('Condition Number After Diagonal: ');
disp(condest(matrix));

clear all;
% Reading matrix after
fileRows = fopen('./MatricesResults/resultAfter2Diag/rows.bin', 'rb');
fileColumns = fopen('./MatricesResults/resultAfter2Diag/columns.bin', 'rb');
fileValues = fopen('./MatricesResults/resultAfter2Diag/values.bin', 'rb');
rows = fread(fileRows, 'int32');
columns = fread(fileColumns, 'int32');
values = fread(fileValues, 'double');
fclose(fileRows);
fclose(fileColumns);
fclose(fileValues);

matrix = sparse(rows, columns, values);

disp('Condition Number After 2 Times Diagonal: ');
disp(condest(matrix));

clear all;
% Reading matrix after
fileRows = fopen('./MatricesResults/resultAfter3Diag/rows.bin', 'rb');
fileColumns = fopen('./MatricesResults/resultAfter3Diag/columns.bin', 'rb');
fileValues = fopen('./MatricesResults/resultAfter3Diag/values.bin', 'rb');
rows = fread(fileRows, 'int32');
columns = fread(fileColumns, 'int32');
values = fread(fileValues, 'double');
fclose(fileRows);
fclose(fileColumns);
fclose(fileValues);

matrix = sparse(rows, columns, values);

disp('Condition Number After 3 Times Diagonal: ');
disp(condest(matrix));

%cond(matrix,1)
% si usas cond (la precisa) con una matriz sparse hace condest por defecto

% Forma precisa de hacer el mismo calculo que cond
%inv_mat = inv(matrix) % Se queda sin memoria al igual que python
%disp(norm(matrix, 1) * norm(inv_mat, 1))