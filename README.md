# Optimized Preconditioners on OpenFOAM

## Resumen

En los últimos años, las simulaciones de dinámica de fluidos (CFD por sus siglas en inglés) han cobrado una gran importancia en la industria. Estas simulaciones realizan costosos cálculos, destacando la resolución de grandes sistemas de ecuaciones mediante el uso de procesos iterativos. Para acelerar dichos procesos es común el uso de unas operaciones intermedias conocidas como Preconditioners, las cuales consiguen restar complejidad gracias a la reducción del número de iteraciones necesarias para converger.

De esta forma, nace un fuerte interés en aumentar el rendimiento de las costosas aplicaciones CFD, tanto en términos de reducción de las horas de cómputo como en términos de consumo energético, motivando así este estudio centrado en la operación conocida como Preconditioner. En concreto, se busca comparar un método sencillo, el Diagonal Preconditioner, frente a otro más complejo pero más beneficioso a priori, el DIC Preconditioner.

Con esto surge un doble objetivo: por un lado, probar empíricamente el beneficio ofrecido por estas operaciones en una aplicación real, OpenFOAM; y, por otro, la paralelización (tanto en CPU como en GPU) de ambos métodos de acondicionamiento, intentando así acelerar la resolución de las simulaciones CFD.

Una vez se concluyeron las pruebas empíricas para determinar la mejora ofrecida por los Preconditioners, se obtuvo un resultado acorde a lo esperado teóricamente. Estas operaciones consiguen acondicionar el sistema actual dotándolo de una mayor capacidad de convergencia, reduciendo así el número de iteraciones necesarias para su resolución.

También se logró la optimización mediante la paralelización, tanto en CPU como en GPU, destacando la gran velocidad conseguida en el método sencillo (Diagonal Preconditioner). Sin embargo, se observa que la menor mejora del método complejo (DIC Preconditioner) puede llegar a ofrecer un rendimiento similar o mejor debido al mayor beneficio visto en el acondicionamiento de la convergencia del sistema a resolver.

## Introducción

La computación de grandes cantidades de datos, así como el cálculo de operaciones matemáticas complejas, han sido siempre un gran reto en el mundo de la informática. A pesar de ello, el esfuerzo por intentar resolver estos problemas es una acción altamente recompensada que puede llevar a grandes avances en áreas como la medicina o la ingeniería.

Hoy en día, el aumento anual de la velocidad de los procesadores ha dejado de ser viable. Como indicaba el cofundador de Intel, Gordon Moore, los procesadores duplicaban el número de transistores cada dos años aproximadamente, impulsando así un crecimiento exponencial en la potencia de cálculo [1]. Sin embargo, esto dejaría de ser una realidad en las últimas décadas por diversas complicaciones físicas, técnicas y económicas, como comenta Esteve Almirall [2]. En este contexto surgen diferentes técnicas para intentar seguir mejorando el rendimiento del cálculo computacional, tales como la optimización de código o recursos, la paralelización de procesos o el uso de arquitecturas específicas [3].

Este trabajo se centra en el uso de técnicas de HPC, en especial la paralelización, para acelerar el proceso de simulación en dinámica de fluidos (CFD por sus siglas en inglés). Nuevamente se enfrentará un problema que requiere de altos costes computacionales, los cuales resultan inviables en sistemas tradicionales, llegando a tardar horas o incluso días en la producción de resultados.

Este tipo de simulaciones trata esencialmente la resolución de grandes sistemas de ecuaciones mediante métodos iterativos (mayor compatibilidad con procesos computacionales). Normalmente, en cada una de las iteraciones que componen este proceso resolutivo, se intercalan unas operaciones conocidas como Preconditioners. Estas operaciones buscan ajustar la matriz que representa al sistema de ecuaciones de tal forma que se reduzca el número de iteraciones necesarias para que dicho proceso converja y, por tanto, muestre una solución en un menor tiempo.

De este modo, el estudio se centra en dos aspectos clave. En primer lugar, evaluar el impacto de los Preconditioners en el sistema a resolver, puesto que hasta la fecha no se han realizado pruebas empíricas al respecto. En segundo lugar, seleccionar dos de estas operaciones, optimizarlas mediante paralelización (tanto en CPU como en GPU) y comparar cuál ofrece mejores resultados en términos de tiempo de cómputo total, considerando tanto el tiempo de ejecución como la reducción en el número de iteraciones.

Para ello, se ha escogido un método rápido y sencillo, el Diagonal Preconditioner, que podría aplicarse varias veces para alcanzar resultados similares a los de un segundo método más complejo, pero con mejor convergencia, el DIC Preconditioner. Además, se han empleado otros cálculos complementarios, como el Condition Number (CN), que permite estimar qué tan cerca está el sistema de ecuaciones de la convergencia.

## Objetivos

Este proyecto busca esencialmente resolver dos hitos. En primer lugar, medir empíricamente la efectividad de la operación conocida como Preconditioner en la resolución de sistemas de ecuaciones mediante métodos iterativos. Para lograrlo, se pretende usar la herramienta matemática conocida como Condition Number, la cual es capaz de resumir lo cerca que se halla la solución actual de la convergencia antes y después de aplicar dicha operación.

En segundo lugar, se busca acelerar la operación mencionada con anterioridad a través de la paralelización. Concretamente, se busca realizar estas pruebas sobre dos tipos de Preconditioners: uno más sencillo, rápido y fácil de paralelizar y otro más complejo y difícil de paralelizar, pero con mejores resultados probados matemáticamente.

Esto permite analizar y comparar distintas soluciones para un mismo problema, teniendo así un enfoque más sencillo, pero que podría ser aplicado varias veces para conseguir los resultados de un sistema más complejo en el caso del Diagonal Preconditioner (opción sencilla). Por otro lado, se analiza un enfoque de mejora de una herramienta compleja, pero con buenos resultados, el Diagonal-based Incomplete Cholesky (DIC) Preconditioner (opción costosa).

De esta forma, la idea principal resulta en el estudio de la mejora que provocan los Preconditioners sobre un sistema de gran tamaño, como los usados en simulaciones CFD. Además, se pretende acelerar y comparar las mejoras proporcionadas por dos tipos de Preconditioners, de forma que se pueda razonar cuál sería la operación más adecuada según el problema que deseemos enfrentar.

## Desarrollo del Proyecto

A continuación se comentan los detalles del desarrollo del proyecto recogido en este repositorio, con objeto de comprender mejor todos aquellos ficheros aquí almacenados que componen el análisis definido en las secciones previas.

### Análisis de Residuos

### Condition Number

### Preconditioners

### Paralelización CPU

### Paralelización GPU

## Resultados 

## Conclusiones
