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
# Parte 3

### Este script ejecuta "Prokka" tres veces con ligeras variaciones para anotar el genoma de un virus. Se especifica el uso de 2 nucleos de CPU, se define un prefijo para los archivos de salida y una etiqueta de locus para las anotaciones. Además, en una de las ejecuciones, se usa una referencia de proteínas para mejorar la anotación, y en otra, se agregan anotaciones adicionales y se forza la sobrescritura de archivos.

Aquí cada paso del script relacionado con la anotación de genes usando **Prokka**.

### Línea 1: 
```bash
prokka --kingdom Viruses --cpus 2 --prefix ERR12389866_prokka --locustag ERR12389866_prokka scaffolds.fasta
```
- **`prokka`**: Ejecuta el programa **Prokka**, que se utiliza para la anotación rápida de genomas.
- **`--kingdom Viruses`**: Especifica que el genoma a anotar pertenece al reino de **virus**. Esto ayuda a Prokka a ajustar los parámetros de anotación para genomas virales.
- **`--cpus 2`**: Utiliza 2 núcleos de CPU para acelerar el proceso de anotación.
- **`--prefix ERR12389866_prokka`**: Especifica el prefijo de los archivos de salida. En este caso, todos los archivos generados tendrán el prefijo `ERR12389866_prokka`.
- **`--locustag ERR12389866_prokka`**: Asigna una etiqueta a los loci anotados (en este caso, `ERR12389866_prokka` será el identificador de los genes).
- **`scaffolds.fasta`**: Es el archivo **FASTA** que contiene las secuencias de andamios (scaffolds) o contigs que se van a anotar.

### Línea 2:
```bash
prokka --kingdom Viruses --cpus 2 --proteins reference.gb --prefix ERR12389866_prokka --locustag ERR12389866_prokka scaffolds.fasta
```
Esta línea es muy similar a la anterior, pero tiene un parámetro adicional:

- **`--proteins reference.gb`**: Indica a Prokka que utilice el archivo **reference.gb** como una referencia de proteínas para la anotación. Este archivo, que debe estar en formato GenBank, contiene proteínas previamente anotadas, que Prokka usará como referencia para mejorar la precisión de la anotación.
- Los demás parámetros son los mismos:
  - **`--kingdom Viruses`**, **`--cpus 2`**, **`--prefix ERR12389866_prokka`**, **`--locustag ERR12389866_prokka`**, y **`scaffolds.fasta`** tienen los mismos significados que en la primera línea.

### Línea 3:
```bash
prokka --kingdom Viruses --cpus 2 --prefix ERR12389866_prokka --locustag ERR12389866_prokka scaffolds.fasta -o ERR12389866_ann --addgenes --force
```
En esta línea se agregan nuevos parámetros adicionales:

- **`-o ERR12389866_ann`**: Define la carpeta de salida para los resultados de la anotación. Aquí, los archivos resultantes se guardarán en el directorio **ERR12389866_ann**.
- **`--addgenes`**: Indica a Prokka que agregue anotaciones de genes en las características anotadas. Esto asegura que cada CDS (secuencia de codificación de proteínas) tenga un gen anotado correspondiente.
- **`--force`**: Sobrescribe cualquier archivo existente en el directorio de salida sin solicitar confirmación.

Los otros parámetros (como **`--kingdom Viruses`**, **`--cpus 2`**, **`--prefix ERR12389866_prokka`**, **`--locustag ERR12389866_prokka`**, y **`scaffolds.fasta`**) son los mismos que en las líneas anteriores.
