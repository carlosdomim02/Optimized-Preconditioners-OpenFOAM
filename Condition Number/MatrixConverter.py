import numpy as np
import scipy as sp
import struct as st

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
        initialResidual = st.unpack(initialResidual_format, file.read(st.calcsize(initialResidual_format)))
        psi = st.unpack(psi_format, file.read(st.calcsize(psi_format)))
        ownerStart = st.unpack(ownerStart_format, file.read(st.calcsize(ownerStart_format)))
    return (diag, lower, lowerAddr, upper, upperAddr, initialResidual, psi, ownerStart)
      

def read_txt_data(file_name):
    with open(file_name) as file:
        diag = list(map(lambda x: float(x), file.readline().split(' ')[:-1]))
        lower = list(map(lambda x: float(x), file.readline().split(' ')[:-1]))
        lowerAddr = list(map(lambda x: int(x), file.readline().split(' ')[:-1]))
        upper = list(map(lambda x: float(x), file.readline().split(' ')[:-1]))
        upperAddr = list(map(lambda x: int(x), file.readline().split(' ')[:-1]))
        initialResidual = list(map(lambda x: float(x), file.readline().split(' ')[:-1]))
        psi = list(map(lambda x: float(x), file.readline().split(' ')[:-1]))
        ownerStart = list(map(lambda x: int(x), file.readline().split(' ')[:-1]))
    return (diag, lower, lowerAddr, upper, upperAddr, initialResidual, psi, ownerStart)


def print_converted_data_bin(result_file, diag, lower, lowerAddr, upper, upperAddr, initialResidual, psi, ownerStart):
    with open(result_file, 'wb') as file:
        file.write(st.pack('8i', len(diag), len(lower), len(lowerAddr), len(upper), len(upperAddr), len(initialResidual), len(psi), len(ownerStart)))
        file.write(st.pack(f'{len(diag)}d', *diag))
        file.write(st.pack(f'{len(lower)}d', *lower))
        file.write(st.pack(f'{len(lowerAddr)}i', *lowerAddr))
        file.write(st.pack(f'{len(upper)}d', *upper))
        file.write(st.pack(f'{len(upperAddr)}i', *upperAddr))
        file.write(st.pack(f'{len(initialResidual)}d', *initialResidual))
        file.write(st.pack(f'{len(psi)}d', *psi))
        file.write(st.pack(f'{len(ownerStart)}i', *ownerStart))


def print_converted_data_txt(result_file, diag, lower, lowerAddr, upper, upperAddr, initialResidual, psi, ownerStart):
    with open(result_file, 'w') as file:
        for value in diag:
            file.write(str(value) + ' ')
        file.write('\n')
        for value in lower:
            file.write(str(value) + ' ')
        file.write('\n')
        for value in lowerAddr:
            file.write(str(value) + ' ')
        file.write('\n')
        for value in upper:
            file.write(str(value) + ' ')
        file.write('\n')
        for value in upperAddr:
            file.write(str(value) + ' ')
        file.write('\n')
        for value in initialResidual:
            file.write(str(value) + ' ')
        file.write('\n')
        for value in psi:
            file.write(str(value) + ' ')
        file.write('\n')
        for value in ownerStart:
            file.write(str(value) + ' ')
        file.write('\n')


def remove_x_vector(file_name, result_file):
    # Read matrix data
    if file_name.endswith('.bin'):
        diag, lower, lowerAddr, upper, upperAddr, initialResidual, psi, ownerStart = read_bin_data(file_name)
    else:
        diag, lower, lowerAddr, upper, upperAddr, initialResidual, psi, ownerStart = read_txt_data(file_name)
    
    # psi (x vector) inverse
    inverse_psi = [1.0/value for value in list(psi)]

    # store for results
    newLower = len(lower)*[0.0]
    newUpper = len(upper)*[0.0]
    
    # A = Ax * psi^-1
    diag = [diag[i] * inverse_psi[i] for i in range(0, len(diag))]
    for k in range(0, len(lowerAddr)):
        i = lowerAddr[k]
        j = upperAddr[k]
        newLower[k] = lower[k] * inverse_psi[i]
        newUpper[k] = upper[k] * inverse_psi[j] 
        #matrix[j,i] = lower[k]
        #matrix[i,j] = upper[k]       
    
    # Printing result
    if result_file.endswith('.bin'):
        print_converted_data_bin(result_file, diag, newLower, lowerAddr, newUpper, upperAddr, initialResidual, psi, ownerStart)
    else:
        print_converted_data_txt(result_file, diag, newLower, lowerAddr, newUpper, upperAddr, initialResidual, psi, ownerStart)
        
    print('matrix proccessing finished')


def add_x_vector(file_name, result_file):
    # Read matrix data
    if file_name.endswith('.bin'):
        diag, lower, lowerAddr, upper, upperAddr, initialResidual, psi, ownerStart = read_bin_data(file_name)
    else:
        diag, lower, lowerAddr, upper, upperAddr, initialResidual, psi, ownerStart = read_txt_data(file_name)

    # store for results
    newLower = len(lower)*[0.0]
    newUpper = len(upper)*[0.0]

    # Ax = A * psi
    diag = [diag[i] * psi[i] for i in range(0, len(diag))]
    for k in range(0, len(lowerAddr)):
        i = lowerAddr[k]
        j = upperAddr[k]
        newLower[k] = lower[k] * psi[i]
        newUpper[k] = upper[k] * psi[j] 
        #matrix[j,i] = lower[k]
        #matrix[i,j] = upper[k]       
    
    # Printing result
    if result_file.endswith('.bin'):
        print_converted_data_bin(result_file, diag, newLower, lowerAddr, newUpper, upperAddr, initialResidual, psi, ownerStart)
    else:
        print_converted_data_txt(result_file, diag, newLower, lowerAddr, newUpper, upperAddr, initialResidual, psi, ownerStart)
        
    print('matrix proccessing finished')


if __name__ == '__main__':
    for itSimpleFoamDir in list(range(100, 501, 100)) + [250]:
        it = str(itSimpleFoamDir)
        remove_x_vector('./MatricesConXIT' + it + '/matrixBefore.bin', './MatricesSinXIT' + it + '/matrixBefore.bin')  
        remove_x_vector('./MatricesConXIT' + it + '/matrixAfterDiag.bin', './MatricesSinXIT' + it + '/matrixAfterDiag.bin')  
        remove_x_vector('./MatricesConXIT' + it + '/matrixAfter2Diag.bin', './MatricesSinXIT' + it + '/matrixAfter2Diag.bin')  
        remove_x_vector('./MatricesConXIT' + it + '/matrixAfter3Diag.bin', './MatricesSinXIT' + it + '/matrixAfter3Diag.bin')  
        remove_x_vector('./MatricesConXIT' + it + '/matrixAfterDIC.bin', './MatricesSinXIT' + it + '/matrixAfterDIC.bin')
    # add_x_vector('./MatricesSinX/matrixBefore.bin', './MatricesConX/resultBefore.bin')
