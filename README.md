# Parte 1

### Este script automatiza el proceso de descarga de archivos SRA, su conversión a archivos FASTQ, la compresión de estos, y la generación de informes de control de calidad usando FastQC.

A continuacion el propósito de cada línea del script paso a paso:

```
1. mkdir fastqc_reports:
   - Crea un directorio llamado `fastqc_reports` donde se guardarán los reportes generados por la herramienta **FastQC**.

2. mkdir sra_files:
   - Crea un directorio llamado `sra_files` para almacenar los archivos **SRA** (Sequence Read Archive) descargados.

3. prefetch --max-size 50G --option-file denv.txt:
   - Utiliza la herramienta **prefetch** para descargar archivos **SRA**. El parámetro `--max-size 50G` establece un tamaño máximo de 50 GB por archivo. La opción `--option-file denv.txt` indica que el archivo `denv.txt` contiene las listas de IDs de los datos **SRA** que deseas descargar.

4. mv */*.sra .:
    - Mueve todos los archivos **SRA** descargados desde los subdirectorios a la carpeta actual. La estructura `*/*.sra` indica que busca los archivos en directorios anidados.

5. fasterq-dump --split-files *.sra :
   - Convierte los archivos **SRA** descargados en archivos **FASTQ**. La opción `--split-files` separa las lecturas en archivos individuales para lecturas pareadas (pair-end).

6. gzip *.fastq:
   - Comprime todos los archivos **FASTQ** generados.

7. fastqc -t 2 *.fastq.gz:
   #Ejecuta la herramienta **FastQC** para realizar el control de calidad en los archivos **FASTQ** comprimidos. La opción `-t 2` indica que se usarán 2 procesadores para ejecutar la tarea de forma más rápida.

8. mv *.sra sra_files/:
   - Mueve los archivos **SRA** descargados a la carpeta **sra_files** que creaste previamente.

9. mv *_fastqc.html *_fastqc.zip fastqc_reports/:
   - Mueve los reportes **HTML** y los archivos comprimidos **ZIP** generados por **FastQC** a la carpeta **fastqc_reports**.

10. ls -lh:
    - Muestra el contenido del directorio actual en formato de lista (`ls -l`) y de una manera legible para humanos (`-h`), mostrando el tamaño de los archivos de forma más comprensible (por ejemplo, en KB, MB, GB).

```
# Parte 2

### Este script descarga e instala las herramientas necesarias para procesar y limpiar lecturas de secuenciación, verifica la calidad de las lecturas procesadas y, finalmente, realiza el ensamblaje del genoma utilizando SPAdes.

A continuacion el propósito de cada línea del script paso a paso:

#### Instalación de programas con "conda":
```
1. conda env -n GenEvol:
   - Crea un nuevo ambiente de conda llamado **GenEvol**. Sin embargo, falta el comando para activarlo o crearlo completamente, ya que normalmente se usaría `conda create -n GenEvol` para crear el ambiente o `conda activate GenEvol` para activarlo.

2. conda install bioconda::fastqc:
   - Instala la herramienta **FastQC** desde el canal **bioconda**, que es un canal especializado en herramientas bioinformáticas. Esto permitirá el control de calidad de las lecturas de secuencias.

3. conda install bioconda::trimmomatic:
   - Instala la herramienta **Trimmomatic**, también desde **bioconda**. **Trimmomatic** se utiliza para limpiar y recortar las lecturas de secuencias, eliminando regiones de baja calidad o adaptadores.
```

#### Limpieza de las lecturas (clean reads):

```
4. trimmomatic PE -threads 2 ERR12389866_1.fastq.gz ERR12389866_2.fastq.gz ERR12389866_1_paired.fq.gz ERR12389866_1_unpaired.fq.gz ERR12389866_2_paired.fq.gz ERR12389866_2_unpaired.fq.gz SLIDINGWINDOW:4:20 MINLEN:50:
   - Ejecuta **Trimmomatic** en modo **PE** (paired-end) para limpiar las lecturas emparejadas:

     ERR12389866_1.fastq.gz: Archivo con las primeras lecturas de las secuencias emparejadas.
     ERR12389866_2.fastq.gz: Archivo con las segundas lecturas de las secuencias emparejadas.
     ERR12389866_1_paired.fq.gz: Salida de lecturas emparejadas limpias para el primer archivo.
     ERR12389866_1_unpaired.fq.gz: Salida de lecturas no emparejadas para el primer archivo.
     ERR12389866_2_paired.fq.gz: Salida de lecturas emparejadas limpias para el segundo archivo.
     ERR12389866_2_unpaired.fq.gz: Salida de lecturas no emparejadas para el segundo archivo.
     SLIDINGWINDOW:4:20: Recorta la lectura cuando la media de calidad en una ventana de 4 bases es inferior a 20.
     MINLEN:50: Descartar lecturas que tengan una longitud inferior a 50 bases después del recorte.
```

#### Control de calidad:
```
5. fastqc -t 2 *_paired.fq.gz:
   - Ejecuta **FastQC** para verificar la calidad de las lecturas emparejadas limpias. El parámetro `-t 2` utiliza dos hilos para procesar más rápidamente.
```

#### Ensamblaje con "SPAdes":
```
6. spades --pe1-1 ERR12389866_1_paired.fq.gz --pe1-2 ERR12389866_2_paired.fq.gz -t 2 -o ERR12389866_2_spades:
   - Ejecuta el programa de ensamblaje **SPAdes** para ensamblar el genoma a partir de las lecturas emparejadas limpias:

     --pe1-1 ERR12389866_1_paired.fq.gz: Archivo de lecturas emparejadas limpias (primer par).
     --pe1-2 ERR12389866_2_paired.fq.gz: Archivo de lecturas emparejadas limpias (segundo par).
     -t 2: Utiliza dos hilos para realizar el ensamblaje.
     -o ERR12389866_2_spades: Define el directorio de salida donde se guardarán los resultados del ensamblaje.

```
