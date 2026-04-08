# Optimized Preconditioners on OpenFOAM

## Resumen

En los últimos años, las simulaciones de dinámica de fluidos (CFD por sus siglas en inglés) han cobrado una gran importancia en la industria. Estas simulaciones realizan costosos cálculos, destacando la resolución de grandes sistemas de ecuaciones mediante el uso de procesos iterativos. Para acelerar dichos procesos es común el uso de unas operaciones intermedias conocidas como Preconditioners, las cuales consiguen restar complejidad gracias a la reducción del número de iteraciones necesarias para converger.

De esta forma, nace un fuerte interés en aumentar el rendimiento de las costosas aplicaciones CFD, tanto en términos de reducción de las horas de cómputo como en términos de consumo energético, motivando así este estudio centrado en la operación conocida como Preconditioner. En concreto, se busca comparar un método sencillo, el Diagonal Preconditioner, frente a otro más complejo pero más beneficioso a priori, el DIC Preconditioner.

Con esto surge un doble objetivo: por un lado, probar empíricamente el beneficio ofrecido por estas operaciones en una aplicación real, OpenFOAM; y, por otro, la paralelización (tanto en CPU como en GPU) de ambos métodos de acondicionamiento, intentando así acelerar la resolución de las simulaciones CFD.

Una vez se concluyeron las pruebas empíricas para determinar la mejora ofrecida por los Preconditioners, se obtuvo un resultado acorde a lo esperado teóricamente. Estas operaciones consiguen acondicionar el sistema actual dotándolo de una mayor capacidad de convergencia, reduciendo así el número de iteraciones necesarias para su resolución.

También se logró la optimización mediante la paralelización, tanto en CPU como en GPU, destacando la gran velocidad conseguida en el método sencillo (Diagonal Preconditioner). Sin embargo, se observa que la menor mejora del método complejo (DIC Preconditioner) puede llegar a ofrecer un rendimiento similar o mejor debido al mayor beneficio visto en el acondicionamiento de la convergencia del sistema a resolver.

## Introducción

La computación de grandes cantidades de datos, así como el cálculo de operaciones matemáticas complejas, han sido siempre un gran reto en el mundo de la informática. A pesar de ello, el esfuerzo por intentar resolver estos problemas es una acción altamente recompensada que puede llevar a grandes avances en áreas como la medicina o la ingeniería.

Hoy en día, el aumento anual de la velocidad de los procesadores ha dejado de ser viable. Como indicaba el cofundador de Intel, Gordon Moore, los procesadores duplicaban el número de transistores cada dos años aproximadamente, impulsando así un crecimiento exponencial en la potencia de cálculo. Sin embargo, esto dejaría de ser una realidad en las últimas décadas por diversas complicaciones físicas, técnicas y económicas, como comenta Esteve Almirall. En este contexto surgen diferentes técnicas para intentar seguir mejorando el rendimiento del cálculo computacional, tales como la optimización de código o recursos, la paralelización de procesos o el uso de arquitecturas específicas.

Este trabajo se centra en el uso de técnicas de HPC, en especial la paralelización, para acelerar el proceso de simulación en dinámica de fluidos (CFD por sus siglas en inglés). Nuevamente se enfrentará un problema que requiere de altos costes computacionales, los cuales resultan inviables en sistemas tradicionales, llegando a tardar horas o incluso días en la producción de resultados.

Este tipo de simulaciones trata esencialmente la resolución de grandes sistemas de ecuaciones mediante métodos iterativos (mayor compatibilidad con procesos computacionales). Normalmente, en cada una de las iteraciones que componen este proceso resolutivo, se intercalan unas operaciones conocidas como Preconditioners. Estas operaciones buscan ajustar la matriz que representa al sistema de ecuaciones de tal forma que se reduzca el número de iteraciones necesarias para que dicho proceso converja y, por tanto, muestre una solución en un menor tiempo.

De este modo, el estudio se centra en dos aspectos clave. En primer lugar, evaluar el impacto de los [`Preconditioners`](https://www.math.iit.edu/~fass/477577_Chapter_16.pdf) en el sistema a resolver, puesto que hasta la fecha no se han realizado pruebas empíricas al respecto. En segundo lugar, seleccionar dos de estas operaciones, optimizarlas mediante paralelización (tanto en CPU como en GPU) y comparar cuál ofrece mejores resultados en términos de tiempo de cómputo total, considerando tanto el tiempo de ejecución como la reducción en el número de iteraciones.

Para ello, se ha escogido un método rápido y sencillo, el [`Diagonal Preconditioner`](https://www.accefyn.com/revista/Vol_28/106/49-55.pdf), que podría aplicarse varias veces para alcanzar resultados similares a los de un segundo método más complejo, pero con mejor convergencia, el [`DIC Preconditioner`](https://docs.nvidia.com/cuda/archive/12.8.0/incomplete-lu-cholesky/index.html). Además, se han empleado otros cálculos complementarios, como el Condition Number (CN), que permite estimar qué tan cerca está el sistema de ecuaciones de la convergencia.

## Objetivos

Este proyecto busca esencialmente resolver dos hitos. En primer lugar, medir empíricamente la efectividad de la operación conocida como Preconditioner en la resolución de sistemas de ecuaciones mediante métodos iterativos. Para lograrlo, se pretende usar la herramienta matemática conocida como Condition Number, la cual es capaz de resumir lo cerca que se halla la solución actual de la convergencia antes y después de aplicar dicha operación.

En segundo lugar, se busca acelerar la operación mencionada con anterioridad a través de la paralelización. Concretamente, se busca realizar estas pruebas sobre dos tipos de Preconditioners: uno más sencillo, rápido y fácil de paralelizar y otro más complejo y difícil de paralelizar, pero con mejores resultados probados matemáticamente.

Esto permite analizar y comparar distintas soluciones para un mismo problema, teniendo así un enfoque más sencillo, pero que podría ser aplicado varias veces para conseguir los resultados de un sistema más complejo en el caso del Diagonal Preconditioner (opción sencilla). Por otro lado, se analiza un enfoque de mejora de una herramienta compleja, pero con buenos resultados, el Diagonal-based Incomplete Cholesky (DIC) Preconditioner (opción costosa).

De esta forma, la idea principal resulta en el estudio de la mejora que provocan los Preconditioners sobre un sistema de gran tamaño, como los usados en simulaciones CFD. Además, se pretende acelerar y comparar las mejoras proporcionadas por dos tipos de Preconditioners, de forma que se pueda razonar cuál sería la operación más adecuada según el problema que deseemos enfrentar.

## Desarrollo del Proyecto

A continuación se comentan los detalles del desarrollo del proyecto recogido en este repositorio, con objeto de comprender mejor todos aquellos ficheros aquí almacenados que componen el análisis definido en las secciones previas. 

Para el análisis realizado a lo largo de este proyecto se escoge un tutorial de OpenFOAM como caso de pruebas, el cual consiste en un sistema de ecuaciones de un caso real, de forma que los resultados obtenidos sean datos realistas de lo que podría suceder en una aplicación real. Concretamente se ha buscado una matriz que represente un sistema de un tamaño considerable en aplicaciones CFD, ya que serán las principales beneficiarias de este estudio. En concreto, en los próximos pasos se usará como matriz de entrada aquella que representa el tutorial “motorBike” de OpenFOAM: [22]

<img width="681" height="422" alt="image" src="https://github.com/user-attachments/assets/eee17c5c-ef5c-4cdf-9e8a-5a9c135f6be1" />

Dicho tutorial es un caso sencillo, pero realista que realiza la simulación de la aerodinámica de una moto. Con esto se obtiene una matriz de dimensiones aproximadas de 300K x 300K elementos. Aunque es posible que haya ejemplos más grandes en problemas reales, este caso resultará lo suficientemente grande para estudiar la eficacia y paralelización de los Preconditioners.

### Análisis de Residuos

En primer lugar, se comienza con la ejecución del caso de pruebas con objeto de medir la supuesta mejora que ofrecen los Preconditioners. Esto permite ver el efecto que refleja su uso en una aplicación real gracias a la información de los “logs” de OpenFOAM. En concreto, se ejecutó la simulación configurando el uso de cada uno de los Preconditioners estudiados. Además, se incluyó una ejecución doble del Diagonal Preconditioner, así como otra ejecución sin ninguna de estas operaciones.

Dentro de estos “logs” proporcionados por OpenFOAM (visibles en [`análisis de residuos`](https://github.com/carlosdomim02/Optimized-Preconditioners-OpenFOAM/tree/main/Residual%20analysis)), destacan medidas como los residuos o el número de iteraciones hasta converger. Los residuos o residual es una medida del desequilibrio o error en las ecuaciones de conservación y se define como la diferencia entre el valor de cierta variable con respecto a esa misma variable obtenida en la iteración anterior [23]. De esta forma, se aplica una técnica muy extendida a la hora de monitorizar este tipo de simulaciones, construir gráficas que permitan analizar dichos residuos u otras medidas de una forma más visual (esto se realiza gracias con la ayuda del script python [`log_analyzer.py`](https://github.com/carlosdomim02/Optimized-Preconditioners-OpenFOAM/blob/main/Residual%20analysis/log_analyzer.py)), las cuales se aprecian en [`SimpleFoamAnalysis.xlsx)`](https://github.com/carlosdomim02/Optimized-Preconditioners-OpenFOAM/blob/main/Residual%20analysis/result%20csvs/SimpleFoamAnalysis.xlsx). 

Con el resumen proporcionado con las gráficas no solo se proporciona información sobre cómo se va acercando el sistema a la convergencia, sino que puede ser de gran ayuda para ver otros detalles como la evolución en cuanto a número de iteraciones necesarias para converger.

### Extracción de los Preconditioners

Es posible que el análisis realizado en el paso anterior este considerando otros cálculos que se hacen dentro de la aplicación, de manera que sería interesante estudiar el problema de forma aislada. Con esto, se extrae todo lo necesario para poder ejecutar los Preconditioners estudiados fuera de OpenFOAM, lo cual se recoge en [`Original Preconditioners`](https://github.com/carlosdomim02/Optimized-Preconditioners-OpenFOAM/tree/main/Original%20Preconditioners).

Esta extracción no solo incluye los algoritmos que realizan las operaciones de acondicionamiento de la matriz, sino que debe obtenerse una o varias matrices que representen un sistema de ecuaciones real. Como se ha mencionado con anterioridad, se toman estos datos del tutorial “motorBike” de OpenFOAM. En concreto, se rescata la matriz que representa la solución actual en varios puntos de la simulación, evitando errores en los análisis debidos a casos puntuales.

En resumen, se extrae todo lo necesario para completar los siguientes pasos fuera de la aplicación OpenFOAM. Desde los algoritmos que deseamos optimizar (el Diagonal Preconditioner y el DIC Preconditioner), hasta los datos de entrada que representan un caso realista.

### Condition Number

En este punto se cuenta con varias matrices que actúan como entrada para los Preconditioners, las cuales pueden ser analizadas desde el punto de vista de la cercanía a la convergencia. Esto representaría el caso de estudio donde aún no se ha aplicado ninguna operación intermedia.

Además, ya se poseen los algoritmos capaces de generar la matriz que representa el sistema tras aplicar cada uno de los Preconditioners sobre los ejemplos originales. De esta manera, es posible conocer la cercanía a la convergencia antes y después de la aplicación de los distintos Preconditioners estudiados.

Por todo ello, es el punto ideal para construir o adaptar una herramienta capaz de calcular el resumen de convergencia ofrecido por el Condition Number sobre los ejemplos objeto de análisis.

Esta herramienta o conjunto de herramientas se centra en el cálculo del [`Condition Number (CN)`](https://www.phys.uconn.edu/~rozman/Courses/m3511_18s/downloads/condnumber.pdf), un método matemático que permite aproximar cómo de cerca se 
encuentra la matriz que representa el sistema actual de converger. El Condition Number, o número de condición en español, tiene como objetivo medir la sensibilidad a pequeñas perturbaciones en los datos de entrada de la matriz, así como a los errores de redondeo que ocurren durante el proceso de solución. Por lo general, un Condition Number no solo se aplica a una matriz en particular, sino también al problema que se está resolviendo. Es posible que una matriz pueda estar mal condicionada para el cálculo de su inversa, mientras que el problema de encontrar sus valores propios está bien condicionado, o viceversa. [18]

En el caso de esta investigación, se usa el CN para comprobar como de bien condicionado está el sistema actual antes y después de aplicar un Preconditioner. En otras palabras, se obtiene un número que resume cómo de cerca está la solución actual de converger antes y después de aplicar la supuesta mejora causada por los Preconditioners estudiados. 

En concreto, se crea el script [`ConditionNumberFrombin.m`](https://github.com/carlosdomim02/Optimized-Preconditioners-OpenFOAM/blob/main/Condition%20Number/ConditionNumberFrombin.m) de MatLab, capaz de calcular el CN para una matriz dada, la cual debe seguir un formato específico que se consigue tras la aplicación de los scripts [`MatrixConverter.py`](https://github.com/carlosdomim02/Optimized-Preconditioners-OpenFOAM/blob/main/Condition%20Number/MatrixConverter.py) y [`MatrixToCondNumFormat.py`](https://github.com/carlosdomim02/Optimized-Preconditioners-OpenFOAM/blob/main/Condition%20Number/MatrixToCondNumFormat.py) de Python (en este mismo orden). Estos scripts buscan separar el estado actual de las variables en el proceso de resolución (OpenFOAM almacena estos valores combinados con la matriz) de las matrices rescatadas de OpenFOAM para realizar las pruebas, así como adaptar estos datos para ser fácilmente legibles por un script MatLab (adaptar el formato de almacenamiento de OpenFOAM a ficheros binarios fácilmente legibles por MatLab).

### Paralelización CPU

El proceso de aceleración de los Preconditioners extraídos de la aplicación OpenFOAM comenzó con el enfoque más sencillo, la paralelización a nivel CPU. Al fin y al cabo, este tipo de paralelización no supone una curva de aprendizaje demasiado agresiva para aquellos que ya conocen la programación y computación clásicas. 

<img width="529" height="369" alt="image" src="https://github.com/user-attachments/assets/1a6ca2b8-dd90-405a-a1f3-10dbe5b4b76b" />

En concreto la aceleración realizada se centra en el uso de multiprocesadores o procesadores multinúcleo de memoria compartida. Estos dispositivos pueden dividir el trabajo en varios procesadores (núcleos), de forma que se elija en cada momento que datos serán compartidos y cuales no a pesar de usar la misma memoria principal todos ellos. Usando también instrucciones vectoriales (AVX), las cuales ejecutan una misma instrucción sobre varios datos de forma simultánea sin ninguna sobrecarga extra.

Para implementar las distintas versiones paralelas de los Preconditioners (estas versiones se pueden visualizar en [`Parallel Preconditioners CPU`](https://github.com/carlosdomim02/Optimized-Preconditioners-OpenFOAM/tree/main/Parallel%20Preconditioners%20CPU)) se ha usado la herramienta OpenMP, ya que ofrece a los programadores de memoria compartida una herramienta para crear aplicaciones paralelas en diversas plataformas. Así como una máquina Ubuntu 22.04.5 LTS sobre un procesador AMD Ryzen Threadripper 2990WX 32 Core Processor con 32 núcleos y hasta 64 hilos como entorno de pruebas.

### Paralelización GPU

Tras concluir el paso anterior se nota que la paralelización es un buen camino a la hora de intentar mejorar el rendimiento de los Preconditioners. De esta forma, se opta por probar los resultados ofrecidos mediante la paralelización con GPU, la cual obtiene grandes mejoras en el rendimiento mediante una paralelización masiva de operaciones sencillas, lo cual encaja a la perfección con la problemática a resolver. 

<img width="846" height="514" alt="image" src="https://github.com/user-attachments/assets/d13b9cd2-88c8-4cc8-b49b-6611c9f5cfe4" />

A diferencia de los procesadores multinúcleo, las unidades de procesamiento gráfico cuentan con miles de núcleos más pequeños que pueden ejecutar tareas sencillas de forma independiente, lo que las hace ideales para cargas de trabajo paralelas. De esta forma, las CPUs se centran en la precisión y el orden, mientras que las GPUs son ideales aplicar operaciones sencillas sobre cantidades masivas de datos en paralelo.

Aunque la función principal de las GPU sigue girando en torno a los gráficos y los elementos visuales cada vez más realistas, las GPUs se han convertido en procesadores paralelos de uso más general capaces de ejecutar muchas tareas de forma simultánea para manejar una variedad cada vez mayor de aplicaciones. 

Es en esta última faceta es donde se sacó provecho de dicho elemento para paralelizar los Preconditioners estudiados, cuyas distintas versiones se observan en [`Parallel Preconditioners GPU`](https://github.com/carlosdomim02/Optimized-Preconditioners-OpenFOAM/tree/main/Parallel%20Preconditioners%20GPU) . Concretamente se ha usado el kit de herramientas CUDA de las tarjetas gráficas NVIDIA. Para poner a pruebas estos algoritmos paralelos se usaron dos entornos de pruebas diferentes:
  - La GPU NVIDIA GeForce GTX 1650 with Max-Q Design del dispositivo Windows local ejecutada a través de Visual Studio 2019 junto al paquete CUDA 12.6 y las opciones de compilación release.
  - La GPU NVIDIA GeForce RTX 3090 integrada en un servidor (contenedor Docker) Ubuntu 22.04.5 LTS administrada mediante PyTorch. Su uso se realizó mediante la compilación con NVCC CUDA 12.1 con la opción -O2 que ofrece características de release.

## Resultados 

### Estudio de la Convergencia

En primer lugar, se visualizan dos gráficas obtenidas a partir de análisis de los residuos obtenidos directamente de la ejecución de los distintos Preconditioners sobre el tutorial “motorBike” de OpenFOAM. Una de ellas compara el número de iteraciones necesarias para converger en cada resolución de matrices (ejecución de PCG) efectuada por el algoritmo de resolución concido como “simpleFoam”:
<img width="1224" height="449" alt="image" src="https://github.com/user-attachments/assets/137d22eb-73ab-4c63-ac66-96367abd0dfb" />

Mientras la otra representa el tiempo de ejecución total (del algoritmo “simpleFoam”) para cada uno de los Preconditioners: 

<img width="1144" height="446" alt="image" src="https://github.com/user-attachments/assets/4835fa74-79a7-4c6d-b04f-73862570a2ee" />

Por otro lado, se aprecia una comparativa de los Condition Numbers obtenidos de forma previa a la ejecución de ningún Preconditioner y tras cada uno de los Preconditioners estudiados (incluyendo una doble y triple ejecución del Diagonal): 

<img width="1441" height="633" alt="image" src="https://github.com/user-attachments/assets/8d0785ad-c712-4f34-9e98-7f9f6ea4748d" />

Esta comparativa se hace por cada una de las matrices extraídas que anteriormente se describieron, las cuales se encuentran en [`Condition Number`](https://github.com/carlosdomim02/Optimized-Preconditioners-OpenFOAM/tree/main/Condition%20Number). 

### Paralelización CPU:

En esta sección se aprecian las mejoras surgidas tras la paralelización en CPU (mediante OMP) de los distintos Preconditioners estudiados. Estos resultados son independientes de la matriz de entrada usada, por tanto, se usa cualquiera de ellas (aunque siempre la misma). Empezando por los resultados de optimizar el Diagonal Preconditioner se obtiene un gráfico que indica cuánto mejora el Diagonal Preconditioner por cantidad de hilos usados respecto de la versión secuencial: 

<img width="701" height="268" alt="image" src="https://github.com/user-attachments/assets/372215a8-9e3d-441a-944c-32255113eaee" />

Por otro lado, se aprecia la siguiente mejora con las diferentes versiones paralelas de DIC Preconditioner (mismo formato, pero con más versiones): 

<img width="853" height="379" alt="image" src="https://github.com/user-attachments/assets/713a4880-4bcd-4dac-a985-a5ff3aacd4ff" />

Por último, se ve un resumen de las mejoras obtenidas con las distintas configuraciones de ambos Preconditioners, donde se muestran únicamente los mejores speedups obtenidos por cada versión de los Preconditioners, junto al número de hilos que dieron dicho tiempo: 

<img width="693" height="171" alt="image" src="https://github.com/user-attachments/assets/31f7ed92-9273-4afc-8543-bea112cf2dcb" />

### Paralelización GPU:

De igual forma que en el apartado previo, se estudian las mejoras obtenidas mediante la paralelización en GPU (a través de CUDA) de los Preconditioners optimizados. En cambio, ahora se analizan dos tipos de mejoras por separado, las cuales se basan en la gráfica utilizada para cada conjunto de pruebas. Nuevamente, se da comienzo por las versiones de Diagonal Preconditioner: 

<img width="1223" height="470" alt="image" src="https://github.com/user-attachments/assets/1db97b82-3f78-4345-9752-08ba1f884560" />

Por otro lado, se muestran los resultados de la paralelización CUDA de las distintas versiones de DIC Preconditioner (mismo formato, pero con más versiones): 

<img width="765" height="937" alt="image" src="https://github.com/user-attachments/assets/7c6eb494-0ace-4ae6-bc4a-7add24966bcc" />

Por último, se ve un resumen de las mejoras obtenidas con las distintas configuraciones de ambos Preconditioners, es decir, se presenta un resumen de los speedups dados por los distintos algoritmos paralelos en CUDA, junto al número de hilos por bloque que produjeron tal resultado: 

<img width="1103" height="308" alt="image" src="https://github.com/user-attachments/assets/2d3caceb-af4d-4dfa-b319-be770e29e8dc" />

## Discusión de Resultados

El análisis de residuos concluye con un resultado que indica claramente el beneficio de usar Preconditioners en los procesos de resolución de sistemas de ecuaciones mediante métodos iterativos. Además, en los resultados relativos al análisis de residuos se aprecia como la repetición del Preconditioner más sencillo de estos es contraproducente, resultando en un número de iteraciones similar al de la versión original sin acondicionar. Esta mejora también se ve reflejada en la reducción de los tiempos de ejecución totales del algoritmo principal, “simpleFoam”. Sin embargo, aquí no se aprecia tanta diferencia entre los dos Preconditioners estudiados como en el número de iteraciones. 

Por otro lado, la teoría que indica que el DIC Preconditioner debería ser un método que condiciona mejor la convergencia de sistema respecto al Diagonal Preconditioner, se ve reforzada por los resultados del estudio del CN. Apreciándose nuevamente la gran mejora ofrecida por los Preconditioners que se aplican una única vez en relación con la versión no acondicionada. También se hace notar el empeoramiento de la convergencia del sistema actual cuando se repite varias veces la aplicación de una de estas operaciones. 

En lo relativo a las mejoras de los Preconditioners optimizados mediante OMP, se aprecia una buena escalabilidad entre todas las versiones, tanto del DIC como el Diagonal Preconditioner. Además, se aprecia una mayor mejora en el método “init” que en “precondition” en la mayoría de los casos. Finalmente, en la tabla resumen se puede ver cómo el método sencillo, Diagonal, resulta en una mayor mejora de rendimiento, mientras que la versión en el orden de las celdas del DIC Preconditioner (DICCells) resulta ser la que consigue una menor mejora. 

Por último, en aquello referente a las optimizaciones mediante CUDA, se aprecia una peor escalabilidad del problema en todas las versiones estudiadas, ya que esta se ve estancada (o incluso se va a un peor rendimiento) a partir de los 32 hilos por bloque independientemente de la gráfica usada para la ejecución de las pruebas. Sin embargo, si se nota una buena escalabilidad de la gráfica menos potente a la más actual, donde los tiempos llegan a ser más de 10 veces mejores. En cuanto a la tabla resumen, se observa como nuevamente el método sencillo es el más beneficiado. En cambio, en este caso, los métodos que respetan el orden de las dependencias (DICUpperCUDA y DICUpperNACUDA) son los que consiguen una menor mejora. 

## Conclusiones

La mejora proporcionada por los Preconditioners, ya no solo queda probada de forma teórica, sino que se consigue probar gracias a este estudio de manera empírica en una aplicación real. Además, no solo ha sido de utilidad para reforzar la teoría, también, ha sido de gran ayuda para ver que el uso repetido de un Preconditioner sencillo (Diagonal) simplemente representa una pérdida de tiempo y no una mayor mejora en la convergencia respecto al método complejo (DIC). Esto hace que la teoría inicial en la que se suponía un posible mejor acondicionamiento en un tiempo de ejecución mejor o similar al repetir dicho método sencillo no sea posible.  

De esta forma, el método DIC refleja una mayor mejora en la convergencia, tal y como se esperaba, aunque el método Diagonal, también consigue una mejora interesante, siendo este método más rápido. Especialmente, esta mayor velocidad de ejecución del método sencillo se ve reflejada en las versiones optimizadas. Era de esperar que este método debido a la falta de dependencias u otros problemas que requiriesen de sincronización, así como la mayor sencillez de sus operaciones, lo hagan no solo más rápido de inicio, sino con mayor capacidad de mejora. 

Por otro lado, en la mejora del método complejo, DIC, se aprecia como las versiones con una sincronización más estricta (DICUpperCUDA y DICUpperNACUDA), las cuales respetan el orden de las dependencias y, por tanto, dan un resultado exactamente igual a las versiones secuenciales; resultan ser más lentas que otras medidas más blandas. De esta forma, se llega a un dilema en el que es posible que la aceleración de la simulación completa sea mayor en el caso más lento, pero estricto, o en el caso más relajado y rápido. Sin embargo, cualquiera de las versiones llegaría a un resultado válido, ya que la aplicación de los Preconditioners consigue una aproximación, la cual puede perder algo de calidad al alterar el orden de las operaciones del DIC, pero no llegaría a dar un resultado erróneo como harían las condiciones de carrera. 
