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

# Parte 4

### Este script automatiza el proceso de descargar datos de NCBI, realizar el mapeo de las secuencias al genoma de referencia, generar archivos BAM, extraer el genoma consenso y anotar los resultados.

A continuacion el propósito de cada línea del script paso a paso:

### 0 Descargar datos de NCBI:
```bash
prefetch --max-size 50G --option-file accessions_mpox.txt
```
- **`prefetch`**: Comando para descargar archivos desde NCBI SRA.
- **`--max-size 50G`**: Limita el tamaño máximo de descarga a 50 GB.
- **`--option-file accessions_mpox.txt`**: Especifica un archivo de texto con los números de acceso SRA que se van a descargar, en este caso para un estudio relacionado con la viruela del mono (mpox).

```bash
mv */*.sra .
```
- Mueve todos los archivos **SRA** descargados de los subdirectorios a la carpeta actual.

```bash
fasterq-dump --split-files *.sra
```
- Convierte los archivos **SRA** a archivos **FASTQ** utilizando **fasterq-dump** y divide los archivos en lecturas **paired-end** (_1.fastq y _2.fastq).

```bash
gzip *fastq
```
- Comprime los archivos **FASTQ** generados en formato **gzip** para ahorrar espacio.

```bash
fastqc *
```
- Ejecuta **FastQC** en todos los archivos **FASTQ** para evaluar la calidad de las lecturas.

---

### 1 Indexar el genoma de referencia:
```bash
bwa index reference.fasta
```
- **`bwa index`**: Indexa el archivo **FASTA** del genoma de referencia, lo que permite una búsqueda rápida de secuencias durante la alineación.
- **`reference.fasta`**: Archivo FASTA que contiene el genoma de referencia al cual se alinearán las lecturas.

---

### 2 Preparar las instrucciones generales:
```bash
for r1 in *fastq.gz
do
prefix=$(basename $r1 _1.fastq.gz)
r2=${prefix}_2.fastq.gz
```
- Este bucle **for** recorre todos los archivos **FASTQ.gz**.
- **`basename $r1 _1.fastq.gz`**: Extrae el prefijo del archivo de lectura **paired-end** (quita `_1.fastq.gz` para obtener el nombre base del archivo).
- **`r2=${prefix}_2.fastq.gz`**: Define la segunda lectura emparejada.

---

### 3 Instrucciones para generar el archivo .bam:
```bash
bwa mem -t 4 reference.fasta $r1 $r2 > ${prefix}_uno.sam
```
- **`bwa mem`**: Alinea las lecturas de **FASTQ** al genoma de referencia utilizando el algoritmo **BWA-MEM**.
- **`-t 4`**: Utiliza 4 núcleos de CPU para acelerar el proceso.
- **`reference.fasta`**: El genoma de referencia.
- **`$r1 $r2`**: Las dos lecturas emparejadas.
- El resultado de la alineación se guarda en un archivo **SAM** llamado `${prefix}_uno.sam`.

```bash
samtools view -@ 4 -bS -T reference.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam
```
- Convierte el archivo **SAM** a **BAM** utilizando **samtools view**. El archivo BAM es un formato comprimido más eficiente.
- **`-@ 4`**: Usa 4 núcleos.
- **`-bS`**: Convierte de **SAM** a **BAM**.
- El archivo de salida es `${prefix}_unoa.bam`.

```bash
samtools sort -@ 4 -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam
```
- Ordena el archivo BAM por nombre de lectura, y lo guarda como `${prefix}_dosa.bam`.

```bash
samtools fixmate -@ 4 -m ${prefix}_dosa.bam ${prefix}_tresa.bam
```
- **`fixmate`** ajusta las lecturas emparejadas en el archivo BAM para que tengan la información correcta.

```bash
samtools sort -@ 4 ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam
```
- Vuelve a ordenar el archivo BAM, esta vez por posición genómica, y lo guarda como `${prefix}_cuatroa.bam`.

```bash
samtools markdup -@ 4 ${prefix}_cuatroa.bam ${prefix}.bam
```
- **`markdup`** marca las lecturas duplicadas en el archivo BAM para poder filtrarlas posteriormente.

```bash
samtools index -@ 4 ${prefix}.bam
```
- Indexa el archivo **BAM** final `${prefix}.bam` para permitir accesos rápidos a las posiciones.

```bash
rm ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ${prefix}_cuatroa.bam
```
- Borra los archivos intermedios para ahorrar espacio.

---

### 4 Extraer genomas consenso:
```bash
for r1 in *bam
do
prefix=$(basename $r1 .bam)
```
- Bucle **for** que recorre todos los archivos **BAM**.

```bash
samtools mpileup -aa -A -d 0 -Q 0 $r1 | ivar consensus -p ${prefix}.fasta -q 25 -t 0.6 -m 10
```
- **`samtools mpileup`** genera un archivo de alineación base por base (mpileup) del archivo BAM.
  - **`-aa`**: Incluye todas las bases.
  - **`-A`**: No ignora las lecturas secundarias.
  - **`-d 0`**: Sin límite en la profundidad de cobertura.
  - **`-Q 0`**: No filtra por calidad de las lecturas.
- **`ivar consensus`**: Genera un genoma consenso basado en la alineación:
  - **`-p ${prefix}.fasta`**: Define el nombre de salida del genoma consenso.
  - **`-q 25`**: Usa lecturas con calidad mayor a 25.
  - **`-t 0.6`**: Para llamar una base, al menos el 60% de las lecturas deben estar de acuerdo.
  - **`-m 10`**: Mínimo de 10 lecturas de cobertura para llamar una base.

---

### 5 Anotación del genoma:
```bash
mkdir -p annotation
mkdir -p ffn
```
- Crea directorios para guardar las anotaciones (**annotation**) y archivos **FFN** (**ffn**).

```bash
for r1 in *fa
do
prefix=$(basename $r1 .fa)
prokka --cpus 4 $r1 -o ${prefix} --prefix ${prefix} --kingdom Viruses
mv ${prefix}/*.gff annotation/${prefix}.gff
done
```
- Bucle **for** que anota cada archivo de genoma en formato **FASTA** (`*.fa`).
- **`prokka`** realiza la anotación:
  - **`--cpus 4`**: Usa 4 hilos de CPU.
  - **`--kingdom Viruses`**: Define el reino como virus.
  - **`-o ${prefix}`**: Define el directorio de salida.
  - **`--prefix ${prefix}`**: Prefijo para los archivos de salida.
- Después de la anotación, mueve el archivo GFF generado a la carpeta **annotation**.

```bash
cp */*.ffn ffn/
```
- Copia los archivos **FFN** (archivos de secuencias de genes anotados en formato FASTA) al directorio **ffn**.

```bash
ls
```
- Lista los archivos en el directorio actual para verificar que todo se haya ejecutado correctamente.

---


# Parte 5

### Mediante el programa Artemis, vamos a identificar las regiones que Prokka ha anotado como CDS (regiones que codifican proteínas). Al explorar estas anotaciones, se puede observar si alguna de las regiones aparece como "hypothetical", lo que podria indicar que no se ha identificado una función conocida para esa región o gen.

A continuacion el propósito de cada línea del script paso a paso:

### Abrir el terminal y usar Artemis:
```bash
conda activate art
art

**Aqui se abrira en una ventana el programa Artemis. Ejecutar:**

Abrir menu file "open file manager" y cargar el archivo del genoma en formato fasta: *fa
luego de cargar el genoma realizar:

Abrir menu file "read and entry" y cargar el archivo en formato gff: *gff

Click derecho en el listado de descripcion de los fragmentos en parte inferior y marcar:
- Show Gene NAmes
- Show Products
Revisar otras caracteristicas
```
- **conda activate art**: Activa el entorno de `conda` llamado `art`, que contiene la herramienta `ARTEMIS` para la visualización y analisis de genomas.
- **art**: Inicia el programa `ARTEMIS`.
- **open file manager\*.fa**: En el menú del programa, selecciona "Archivo" y luego "Abrir el administrador de archivos" para cargar el archivo del genoma en formato FASTA (`.fa`).
- **read and entry \*.gff**: En el menú del programa, selecciona "Archivo" y luego "Leer y entrada" para cargar el archivo de anotación GFF, que contiene información sobre las regiones genómicas anotadas, como genes y CDS (códigos para proteínas).
- **Show Gene Names**: Para visualizar los nombres de los genes identificados.
- **Show Products**: Para mostrar los productos génicos (como proteínas).

![Descripción de la imagen](https://github.com/Vidal-Felices/Gen_Evol_VFA/raw/main/Captura%20de%20pantalla%20de%202024-10-07%2010-43-15.png)
