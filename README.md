# Optimized Preconditioners on OpenFOAM

## Introducción

En los últimos años, las simulaciones de dinámica de fluidos (CFD por sus siglas en inglés) han cobrado una gran importancia en la industria. Estas simulaciones realizan costosos cálculos, destacando la resolución de grandes sistemas de ecuaciones mediante el uso de procesos iterativos. Para acelerar dichos procesos es común el uso de unas operaciones intermedias conocidas como Preconditioners, las cuales consiguen restar complejidad gracias a la reducción del número de iteraciones necesarias para converger.

De esta forma, nace un fuerte interés en aumentar el rendimiento de las costosas aplicaciones CFD, tanto en términos de reducción de las horas de cómputo como en términos de consumo energético, motivando así este estudio centrado en la operación conocida como Preconditioner. En concreto, se busca comparar un método sencillo, el Diagonal Preconditioner, frente a otro más complejo pero más beneficioso a priori, el DIC Preconditioner.

Con esto surge un doble objetivo: por un lado, probar empíricamente el beneficio ofrecido por estas operaciones en una aplicación real, OpenFOAM; y, por otro, la paralelización (tanto en CPU como en GPU) de ambos métodos de acondicionamiento, intentando así acelerar la resolución de las simulaciones CFD.

Una vez se concluyeron las pruebas empíricas para determinar la mejora ofrecida por los Preconditioners, se obtuvo un resultado acorde a lo esperado teóricamente. Estas operaciones consiguen acondicionar el sistema actual dotándolo de una mayor capacidad de convergencia, reduciendo así el número de iteraciones necesarias para su resolución.

También se logró la optimización mediante la paralelización, tanto en CPU como en GPU, destacando la gran velocidad conseguida en el método sencillo (Diagonal Preconditioner). Sin embargo, se observa que la menor mejora del método complejo (DIC Preconditioner) puede llegar a ofrecer un rendimiento similar o mejor debido al mayor beneficio visto en el acondicionamiento de la convergencia del sistema a resolver.

## Objetivos


## Desarrollo del Proyecto

### Análisis de Residuos

### Condition Number

### Preconditioners

### Paralelización CPU

### Paralelización GPU

## Resultados 

## Conclusiones
