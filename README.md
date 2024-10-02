# Parte 1

### Este script automatiza el proceso de descarga de archivos SRA, su conversión a archivos FASTQ, la compresión de estos, y la generación de informes de control de calidad usando FastQC.

A continuacion el propósito de cada línea de tu script paso a paso:

```
1. mkdir fastqc_reports:
   #Crea un directorio llamado `fastqc_reports` donde se guardarán los reportes generados por la herramienta **FastQC**.

2. mkdir sra_files:
   #Crea un directorio llamado `sra_files` para almacenar los archivos **SRA** (Sequence Read Archive) descargados.

3. prefetch --max-size 50G --option-file denv.txt:
   #Utiliza la herramienta **prefetch** para descargar archivos **SRA**. El parámetro `--max-size 50G` establece un tamaño máximo de 50 GB por archivo. La opción `--option-file denv.txt` indica que el archivo `denv.txt` contiene las listas de IDs de los datos **SRA** que deseas descargar.

4. mv */*.sra .:
   #Mueve todos los archivos **SRA** descargados desde los subdirectorios a la carpeta actual. La estructura `*/*.sra` indica que busca los archivos en directorios anidados.

5. fasterq-dump --split-files *.sra :
   #Convierte los archivos **SRA** descargados en archivos **FASTQ**. La opción `--split-files` separa las lecturas en archivos individuales para lecturas pareadas (pair-end).

6. gzip *.fastq:
   #Comprime todos los archivos **FASTQ** generados.

7. fastqc -t 2 *.fastq.gz:
   #Ejecuta la herramienta **FastQC** para realizar el control de calidad en los archivos **FASTQ** comprimidos. La opción `-t 2` indica que se usarán 2 procesadores para ejecutar la tarea de forma más rápida.

8. mv *.sra sra_files/:
   #Mueve los archivos **SRA** descargados a la carpeta **sra_files** que creaste previamente.

9. mv *_fastqc.html *_fastqc.zip fastqc_reports/:
   #Mueve los reportes **HTML** y los archivos comprimidos **ZIP** generados por **FastQC** a la carpeta **fastqc_reports**.

10. ls -lh:
    #Muestra el contenido del directorio actual en formato de lista (`ls -l`) y de una manera legible para humanos (`-h`), mostrando el tamaño de los archivos de forma más comprensible (por ejemplo, en KB, MB, GB).

```
