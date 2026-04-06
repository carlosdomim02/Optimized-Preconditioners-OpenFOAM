import numpy as np
import scipy as sp
import struct as st

#class MatrixWrapper:
#    def __init__(self, sp_matrix):
#        self.sp_matrix = sp_matrix
#    
#    def __getitem__(self, index):
#        return self.sp_matrix.getrow(index).toarray()[0]


def read_bin_data(file_name):
    with open(file_name, 'rb') as file:
        diag_format = f'{st.unpack("1i", file.read(st.calcsize("1i")))[0]}d'
        lower_format = f'{st.unpack("1i", file.read(st.calcsize("1i")))[0]}d'
        lowerAddr_format = f'{st.unpack("1i", file.read(st.calcsize("1i")))[0]}i'
        upper_format = f'{st.unpack("1i", file.read(st.calcsize("1i")))[0]}d'
        upperAddr_format = f'{st.unpack("1i", file.read(st.calcsize("1i")))[0]}i'
        initialResidual_format = f'{st.unpack("1i", file.read(st.calcsize("1i")))[0]}d'
        psi_format = f'{st.unpack("1i", file.read(st.calcsize("1i")))[0]}d'
        ownerStart_format = f'{st.unpack("1i", file.read(st.calcsize("1i")))[0]}i'

        diag = st.unpack(diag_format, file.read(st.calcsize(diag_format)))
        lower = st.unpack(lower_format, file.read(st.calcsize(lower_format)))
        lowerAddr = st.unpack(lowerAddr_format, file.read(st.calcsize(lowerAddr_format)))
        upper = st.unpack(upper_format, file.read(st.calcsize(upper_format)))
        upperAddr = st.unpack(upperAddr_format, file.read(st.calcsize(upperAddr_format)))
    return (diag, lower, lowerAddr, upper, upperAddr)
      

def read_txt_data(file_name):
    with open(file_name) as file:
        diag = list(map(lambda x: float(x), file.readline().split(' ')[:-1]))
        lower = list(map(lambda x: float(x), file.readline().split(' ')[:-1]))
        lowerAddr = list(map(lambda x: int(x), file.readline().split(' ')[:-1]))
        upper = list(map(lambda x: float(x), file.readline().split(' ')[:-1]))
        upperAddr = list(map(lambda x: int(x), file.readline().split(' ')[:-1]))
    return (diag, lower, lowerAddr, upper, upperAddr)


def print_converted_data_bin(result_folder_name, rows, columns, values):
    with open(f'{result_folder_name}/rows.bin', 'wb') as file:
        file.write(st.pack(f'{len(rows)}i', *rows))

    with open(f'{result_folder_name}/columns.bin', 'wb') as file:
        file.write(st.pack(f'{len(columns)}i', *columns))

    with open(f'{result_folder_name}/values.bin', 'wb') as file:
        file.write(st.pack(f'{len(values)}d', *values))


def print_converted_data_txt(result_folder_name, rows, columns, values):
    with open(f'{result_folder_name}/rows.txt', 'w') as file:
        for row in rows:
            file.write(str(row) + ' ')

    with open(f'{result_folder_name}/columns.txt', 'w') as file:
        for column in columns:
            file.write(str(column) + ' ') 

    with open(f'{result_folder_name}/values.txt', 'w') as file:
        for value in values:
            file.write(str(value) + ' ')


def condition_number_from_file(file_name, result_folder_name, result_type):
    # Read matrix data
    if file_name.endswith('.bin'):
        diag, lower, lowerAddr, upper, upperAddr = read_bin_data(file_name)
    else:
        diag, lower, lowerAddr, upper, upperAddr = read_txt_data(file_name)
    
    # Create and fill matrix
    rows = [k + 1 for k in range(0,len(diag))]
    columns = [k + 1 for k in range(0,len(diag))]
    values = [value for value in diag]

    for k in range(0, len(lowerAddr)):
        i = lowerAddr[k]
        j = upperAddr[k]
        rows += [j + 1]
        columns += [i + 1]
        values += [lower[k]]
        rows += [i + 1]
        columns += [j + 1]
        values += [upper[k]]
        #matrix[j,i] = lower[k]
        #matrix[i,j] = upper[k]            
    
    # Printing result
    if result_type.endswith('bin'):
        print_converted_data_bin(result_folder_name, rows, columns, values)
    else:
        print_converted_data_txt(result_folder_name, rows, columns, values)
        
    print('matrix proccessing finished')


if __name__ == '__main__':
    # Only 1 iteration for CN computation in Matlab
    it = '250'
    condition_number_from_file('./MatricesSinXIT' + it + '/matrixBefore.bin', './MatricesResults/resultBefore', '.bin')    
    condition_number_from_file('./MatricesSinXIT' + it + '/matrixAfterDiag.bin', './MatricesResults/resultAfterDiag', '.bin')
    condition_number_from_file('./MatricesSinXIT' + it + '/matrixAfter2Diag.bin', './MatricesResults/resultAfter2Diag', '.bin')
    condition_number_from_file('./MatricesSinXIT' + it + '/matrixAfter3Diag.bin', './MatricesResults/resultAfter3Diag', '.bin')
    condition_number_from_file('./MatricesSinXIT' + it + '/matrixAfterDIC.bin', './MatricesResults/resultAfterDIC', '.bin')
    