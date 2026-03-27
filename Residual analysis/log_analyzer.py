def analizador_csv(file_name, output_file_name):
    with open(file_name) as file:
        it = 1
        out_csv = ['Iteration,Time,Initial residual,Final residual,N_iterations,Exc_time\n']
        # Select important data
        important_data = [line for line in file 
                    if line.startswith('DICPCG:') or line.startswith('nonePCG:')  or line.startswith("diagonalPCG") or line.startswith('doubleDiagonalPCG:')
                    or line.startswith('Time = ') or line.startswith('ExecutionTime = ')]
        # Group data
        important_data = [tuple(important_data[i:i+3]) for i in range(0,len(important_data),3)]
        # Procesing data
        for time, residual, exc_time in important_data:
            initial, final, its = residual.split(',')[1:]
            out_csv += [(f"{it},"
                         f"{time[time.find('=')+2:time.find('s')]},"
                         f"{initial[initial.find('=')+2:]},"
                         f"{final[final.find('=')+2:]},"
                         f"{its[its.find('s')+2:-1]},"
                         f"{exc_time[exc_time.find('=')+2:exc_time.find('s')]}\n")]
            it += 1
    
    with open(output_file_name, 'w') as file:
        # Write data as csv
        file.writelines(out_csv)

if __name__ == '__main__':
    analizador_csv('./original logs/log.simpleFoam', './result csvs/simpleFoam.csv')
    analizador_csv('./original logs/log.simpleFoamDoubleDiagonal', './result csvs/simpleFoamDoubleDiagonal.csv')
    analizador_csv('./original logs/log.simpleFoamDiagonal', './result csvs/simpleFoamDiagonal.csv')
    analizador_csv('./original logs/log.simpleFoamDic', './result csvs/simpleFoamDic.csv')
    