# Gen_Evol_VFA
A script collection for Evolutive genomics
<pre>(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/SEMESTRE_2024II/G_Evolutiva/Clase2</b></font>$ ls -lh
total 36K
-rw-rw-r-- 1 vidal vidal 1,3K set  9 21:06  clase2.txt
-rw-rw-r-- 1 vidal vidal 1,3K set  7 16:55  commands_1.sh
-rw-rw-r-- 1 vidal vidal 1,3K set  7 22:13  commands_2.sh
-rw-rw-r-- 1 vidal vidal 6,8K set  9 20:53 &apos;Documento sin título.docx&apos;
drwxrwxr-x 2 vidal vidal 4,0K set 15 18:43  <font color="#729FCF"><b>ejemplo_1</b></font>
drwxrwxr-x 2 vidal vidal 4,0K set 15 18:43  <font color="#729FCF"><b>ejemplo_2</b></font>
-rw-rw-r-- 1 vidal vidal   24 set  9 15:49  sra_accessions_1.txt
-rw-rw-r-- 1 vidal vidal   24 set  7 21:29  sra_accessions_2.txt
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/SEMESTRE_2024II/G_Evolutiva/Clase2</b></font>$ ^C
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/SEMESTRE_2024II/G_Evolutiva/Clase2</b></font>$ prefetch --max-size 50G --option-file sra_accessions_1.txt ;
No se ha encontrado la orden «prefetch», pero se puede instalar con:
sudo apt install sra-toolkit
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/SEMESTRE_2024II/G_Evolutiva/Clase2</b></font>$ cd ..
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/SEMESTRE_2024II/G_Evolutiva</b></font>$ cd
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ ls
<font color="#729FCF"><b>anaconda3</b></font>  <font color="#729FCF"><b>Descargas</b></font>  <font color="#729FCF"><b>Documentos</b></font>  <font color="#729FCF"><b>Escritorio</b></font>  <font color="#729FCF"><b>Imágenes</b></font>  <font color="#729FCF"><b>Música</b></font>  <font color="#729FCF"><b>Plantillas</b></font>  <font color="#729FCF"><b>Público</b></font>  <font color="#729FCF"><b>sratoolkit.3.1.1-ubuntu64</b></font>  <font color="#8AE234"><b>stk.tar.gz</b></font>  <font color="#729FCF"><b>Vídeos</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ cd sratoolkit.3.1.1-ubuntu64/
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/sratoolkit.3.1.1-ubuntu64</b></font>$ ls
<font color="#729FCF"><b>bin</b></font>  CHANGES  <font color="#729FCF"><b>example</b></font>  README-blastn  README.md  README-vdb-config  <font color="#729FCF"><b>schema</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/sratoolkit.3.1.1-ubuntu64</b></font>$ prefetch --max-size 50G --option-file sra_accessions_1.txt ;
No se ha encontrado la orden «prefetch», pero se puede instalar con:
sudo apt install sra-toolkit
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/sratoolkit.3.1.1-ubuntu64</b></font>$ conda activate prefetch

EnvironmentNameNotFound: Could not find conda environment: prefetch
You can list all discoverable environments with `conda info --envs`.


(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/sratoolkit.3.1.1-ubuntu64</b></font>$ conda info --envs
# conda environments:
#
base                  *  /home/vidal/anaconda3

(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/sratoolkit.3.1.1-ubuntu64</b></font>$ cd ..
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ export PATH=$PATH:$PWD/sratoolkit.3.0.10-ubuntu64/bin
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ ls -lh
total 90M
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ prefetch --max-size 50G --option-file sra_accessions_1.txt ;
No se ha encontrado la orden «prefetch», pero se puede instalar con:
sudo apt install sra-toolkit
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ ^C
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ prefetch --max-size 50G --option-file sra_accessions_1.txt ;
No se ha encontrado la orden «prefetch», pero se puede instalar con:
sudo apt install sra-toolkit
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ prefetch --max-size 50G --option-file sra_accessions_1.txt ;
No se ha encontrado la orden «prefetch», pero se puede instalar con:
sudo apt install sra-toolkit
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ sudo apt install sra-toolkit
[sudo] contraseña para vidal:            
Leyendo lista de paquetes... Hecho
Creando árbol de dependencias... Hecho
Leyendo la información de estado... Hecho
Se instalarán los siguientes paquetes adicionales:
  blends-common libaec0 libhdf5-103-1t64 libmbedtls14t64 libmbedx509-1t64 libncbi-vdb3 libncbi-wvdb3 libre2-10 libsz2 med-config menu ncbi-vdb-data
Paquetes sugeridos:
  blends-doc menu-l10n gksu | kde-cli-tools | ktsuss
Se instalarán los siguientes paquetes NUEVOS:
  blends-common libaec0 libhdf5-103-1t64 libmbedtls14t64 libmbedx509-1t64 libncbi-vdb3 libncbi-wvdb3 libre2-10 libsz2 med-config menu ncbi-vdb-data
  sra-toolkit
0 actualizados, 13 nuevos se instalarán, 0 para eliminar y 0 no actualizados.
Se necesita descargar 11,3 MB de archivos.
Se utilizarán 30,7 MB de espacio de disco adicional después de esta operación.
¿Desea continuar? [S/n] s
Des:1 http://archive.ubuntu.com/ubuntu noble/universe amd64 menu amd64 2.1.50 [377 kB]
Des:2 http://archive.ubuntu.com/ubuntu noble/universe amd64 blends-common all 0.7.6 [15,2 kB]
Des:3 http://archive.ubuntu.com/ubuntu noble/universe amd64 libaec0 amd64 1.1.2-1build1 [22,9 kB]
Des:4 http://archive.ubuntu.com/ubuntu noble/universe amd64 libsz2 amd64 1.1.2-1build1 [5.476 B]
Des:5 http://archive.ubuntu.com/ubuntu noble/universe amd64 libhdf5-103-1t64 amd64 1.10.10+repack-3.1ubuntu4 [1.270 kB]
Des:6 http://archive.ubuntu.com/ubuntu noble/universe amd64 libmbedx509-1t64 amd64 2.28.8-1 [46,6 kB]
Des:7 http://archive.ubuntu.com/ubuntu noble/universe amd64 libmbedtls14t64 amd64 2.28.8-1 [82,2 kB]
Des:8 http://archive.ubuntu.com/ubuntu noble/universe amd64 ncbi-vdb-data all 3.0.2+dfsg-2build2 [67,5 kB]
Des:9 http://archive.ubuntu.com/ubuntu noble/universe amd64 libncbi-vdb3 amd64 3.0.2+dfsg-2build2 [1.112 kB]
Des:10 http://archive.ubuntu.com/ubuntu noble/universe amd64 libncbi-wvdb3 amd64 3.0.2+dfsg-2build2 [1.137 kB]
Des:11 http://archive.ubuntu.com/ubuntu noble/main amd64 libre2-10 amd64 20230301-3build1 [168 kB]
Des:12 http://archive.ubuntu.com/ubuntu noble/universe amd64 sra-toolkit amd64 3.0.3+dfsg-6ubuntu5 [7.022 kB]
Des:13 http://archive.ubuntu.com/ubuntu noble/universe amd64 med-config all 3.8.1 [11,2 kB]
Descargados 11,3 MB en 5s (2.205 kB/s)<font color="#C4A000"> </font>
Preconfigurando paquetes ...
Seleccionando el paquete menu previamente no seleccionado.
(Leyendo la base de datos ... 497792 ficheros o directorios instalados actualmente.)
Preparando para desempaquetar .../00-menu_2.1.50_amd64.deb ...
Desempaquetando menu (2.1.50) ...
Seleccionando el paquete blends-common previamente no seleccionado.
Preparando para desempaquetar .../01-blends-common_0.7.6_all.deb ...
Desempaquetando blends-common (0.7.6) ...
Seleccionando el paquete libaec0:amd64 previamente no seleccionado.
Preparando para desempaquetar .../02-libaec0_1.1.2-1build1_amd64.deb ...
Desempaquetando libaec0:amd64 (1.1.2-1build1) ...
Seleccionando el paquete libsz2:amd64 previamente no seleccionado.
Preparando para desempaquetar .../03-libsz2_1.1.2-1build1_amd64.deb ...
Desempaquetando libsz2:amd64 (1.1.2-1build1) ...
Seleccionando el paquete libhdf5-103-1t64:amd64 previamente no seleccionado.
Preparando para desempaquetar .../04-libhdf5-103-1t64_1.10.10+repack-3.1ubuntu4_amd64.deb ...
Desempaquetando libhdf5-103-1t64:amd64 (1.10.10+repack-3.1ubuntu4) ...
Seleccionando el paquete libmbedx509-1t64:amd64 previamente no seleccionado.
Preparando para desempaquetar .../05-libmbedx509-1t64_2.28.8-1_amd64.deb ...
Desempaquetando libmbedx509-1t64:amd64 (2.28.8-1) ...
Seleccionando el paquete libmbedtls14t64:amd64 previamente no seleccionado.
Preparando para desempaquetar .../06-libmbedtls14t64_2.28.8-1_amd64.deb ...
Desempaquetando libmbedtls14t64:amd64 (2.28.8-1) ...
Seleccionando el paquete ncbi-vdb-data previamente no seleccionado.
Preparando para desempaquetar .../07-ncbi-vdb-data_3.0.2+dfsg-2build2_all.deb ...
Desempaquetando ncbi-vdb-data (3.0.2+dfsg-2build2) ...
Seleccionando el paquete libncbi-vdb3:amd64 previamente no seleccionado.
Preparando para desempaquetar .../08-libncbi-vdb3_3.0.2+dfsg-2build2_amd64.deb ...
Desempaquetando libncbi-vdb3:amd64 (3.0.2+dfsg-2build2) ...
Seleccionando el paquete libncbi-wvdb3:amd64 previamente no seleccionado.
Preparando para desempaquetar .../09-libncbi-wvdb3_3.0.2+dfsg-2build2_amd64.deb ...
Desempaquetando libncbi-wvdb3:amd64 (3.0.2+dfsg-2build2) ...
Seleccionando el paquete libre2-10:amd64 previamente no seleccionado.
Preparando para desempaquetar .../10-libre2-10_20230301-3build1_amd64.deb ...
Desempaquetando libre2-10:amd64 (20230301-3build1) ...
Seleccionando el paquete sra-toolkit previamente no seleccionado.
Preparando para desempaquetar .../11-sra-toolkit_3.0.3+dfsg-6ubuntu5_amd64.deb ...
Desempaquetando sra-toolkit (3.0.3+dfsg-6ubuntu5) ...
Seleccionando el paquete med-config previamente no seleccionado.
Preparando para desempaquetar .../12-med-config_3.8.1_all.deb ...
Desempaquetando med-config (3.8.1) ...
Configurando libmbedx509-1t64:amd64 (2.28.8-1) ...
Configurando libre2-10:amd64 (20230301-3build1) ...
Configurando libaec0:amd64 (1.1.2-1build1) ...
Configurando libmbedtls14t64:amd64 (2.28.8-1) ...
Configurando ncbi-vdb-data (3.0.2+dfsg-2build2) ...
Configurando menu (2.1.50) ...
Configurando libsz2:amd64 (1.1.2-1build1) ...
Configurando libncbi-vdb3:amd64 (3.0.2+dfsg-2build2) ...
Configurando libncbi-wvdb3:amd64 (3.0.2+dfsg-2build2) ...
Configurando libhdf5-103-1t64:amd64 (1.10.10+repack-3.1ubuntu4) ...
Configurando sra-toolkit (3.0.3+dfsg-6ubuntu5) ...
Procesando disparadores para libc-bin (2.39-0ubuntu8.3) ...
Procesando disparadores para man-db (2.12.0-4build2) ...
Procesando disparadores para install-info (7.1-3build2) ...
Procesando disparadores para menu (2.1.50) ...
Configurando blends-common (0.7.6) ...
Configurando med-config (3.8.1) ...
info: Seleccionando un GID del rango 100 a 999 ...
info: Añadiendo el grupo `med&apos; (GID 127) ...
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ prefetch --max-size 50G --option-file sra_accessions_1.txt ;
2024-09-16T00:03:33 prefetch.3.0.3 int: file not found while opening manager within virtual file system module - VFSManagerOpenFileRead() failed
Usage:
  prefetch [options] &lt;SRA accession&gt; [...]
  Download SRA files and their dependencies

  prefetch [options] --cart &lt;kart file&gt;
  Download cart file

  prefetch [options] &lt;URL&gt; --output-file &lt;FILE&gt;
  Download URL to FILE

  prefetch [options] &lt;URL&gt; [...] --output-directory &lt;DIRECTORY&gt;
  Download URL or URL-s to DIRECTORY

  prefetch [options] &lt;SRA file&gt; [...]
  Check SRA file for missed dependencies and download them

(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ ls
<font color="#729FCF"><b>anaconda3</b></font>  <font color="#729FCF"><b>Descargas</b></font>  <font color="#729FCF"><b>Documentos</b></font>  <font color="#729FCF"><b>Escritorio</b></font>  <font color="#729FCF"><b>Imágenes</b></font>  <font color="#729FCF"><b>Música</b></font>  <font color="#729FCF"><b>Plantillas</b></font>  <font color="#729FCF"><b>Público</b></font>  <font color="#729FCF"><b>sratoolkit.3.1.1-ubuntu64</b></font>  <font color="#8AE234"><b>stk.tar.gz</b></font>  <font color="#729FCF"><b>Vídeos</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ cd Documentos
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos</b></font>$ ls
<font color="#729FCF"><b>UNMSM</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos</b></font>$ cd UNMSM
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM</b></font>$ LS
LS: no se encontró la orden
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM</b></font>$ ls
<font color="#729FCF"><b>SEMESTRE_2024II</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM</b></font>$ cd 2024ii
bash: cd: 2024ii: No existe el archivo o el directorio
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM</b></font>$ cd 2024II
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II</b></font>$ cd G_evol
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol</b></font>$ cd Clase2
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol/Clase2</b></font>$ prefetch --max-size 50G --option-file sra_accessions_1.txt ;

2024-09-16T00:06:03 prefetch.3.0.3: Current preference is set to retrieve SRA Normalized Format files with full base quality scores.
2024-09-16T00:06:03 prefetch.3.0.3: 1) Downloading &apos;ERR12389866&apos;...
2024-09-16T00:06:03 prefetch.3.0.3: SRA Normalized Format file is being retrieved, if this is different from your preference, it may be due to current file availability.
2024-09-16T00:06:03 prefetch.3.0.3:  Downloading via HTTPS...
2024-09-16T00:06:05 prefetch.3.0.3:  HTTPS download succeed
2024-09-16T00:06:05 prefetch.3.0.3:  &apos;ERR12389866&apos; is valid
2024-09-16T00:06:05 prefetch.3.0.3: 1) &apos;ERR12389866&apos; was downloaded successfully
2024-09-16T00:06:05 prefetch.3.0.3: &apos;ERR12389866&apos; has 0 unresolved dependencies

2024-09-16T00:06:06 prefetch.3.0.3: Current preference is set to retrieve SRA Normalized Format files with full base quality scores.
2024-09-16T00:06:06 prefetch.3.0.3: 2) Downloading &apos;ERR12543675&apos;...
2024-09-16T00:06:06 prefetch.3.0.3: SRA Normalized Format file is being retrieved, if this is different from your preference, it may be due to current file availability.
2024-09-16T00:06:06 prefetch.3.0.3:  Downloading via HTTPS...
2024-09-16T00:06:23 prefetch.3.0.3:  HTTPS download succeed
2024-09-16T00:06:23 prefetch.3.0.3:  &apos;ERR12543675&apos; is valid
2024-09-16T00:06:23 prefetch.3.0.3: 2) &apos;ERR12543675&apos; was downloaded successfully
2024-09-16T00:06:23 prefetch.3.0.3: &apos;ERR12543675&apos; has 0 unresolved dependencies
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol/Clase2</b></font>$ cd
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ ls
<font color="#729FCF"><b>anaconda3</b></font>  <font color="#729FCF"><b>Descargas</b></font>  <font color="#729FCF"><b>Documentos</b></font>  <font color="#729FCF"><b>Escritorio</b></font>  <font color="#729FCF"><b>Imágenes</b></font>  <font color="#729FCF"><b>Música</b></font>  <font color="#729FCF"><b>Plantillas</b></font>  <font color="#729FCF"><b>Público</b></font>  <font color="#729FCF"><b>sratoolkit.3.1.1-ubuntu64</b></font>  <font color="#8AE234"><b>stk.tar.gz</b></font>  <font color="#729FCF"><b>Vídeos</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ cd sratoolkit.3.1.1-ubuntu64/
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/sratoolkit.3.1.1-ubuntu64</b></font>$ ls
<font color="#729FCF"><b>bin</b></font>  CHANGES  <font color="#729FCF"><b>example</b></font>  README-blastn  README.md  README-vdb-config  <font color="#729FCF"><b>schema</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/sratoolkit.3.1.1-ubuntu64</b></font>$ cd example
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/sratoolkit.3.1.1-ubuntu64/example</b></font>$ ls
<font color="#729FCF"><b>perl</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/sratoolkit.3.1.1-ubuntu64/example</b></font>$ cd perl
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/sratoolkit.3.1.1-ubuntu64/example/perl</b></font>$ ls
<font color="#8AE234"><b>base-stats.pl</b></font>  <font color="#8AE234"><b>dump-reference.pl</b></font>  <font color="#8AE234"><b>gene-lookup.pl</b></font>  <font color="#8AE234"><b>mismatch-stats.pl</b></font>  <font color="#8AE234"><b>quality-stats.pl</b></font>  <font color="#8AE234"><b>simplefastq.pl</b></font>  <font color="#8AE234"><b>splitfastq.pl</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/sratoolkit.3.1.1-ubuntu64/example/perl</b></font>$ cd ..
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/sratoolkit.3.1.1-ubuntu64/example</b></font>$ cd ..
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/sratoolkit.3.1.1-ubuntu64</b></font>$ ls -lh
total 72K
drwxrwxr-x 3 vidal vidal 4,0K may 13 21:32 <font color="#729FCF"><b>bin</b></font>
-rw-rw-r-- 1 vidal vidal  35K may 13 21:25 CHANGES
drwxrwxr-x 3 vidal vidal 4,0K may 13 21:25 <font color="#729FCF"><b>example</b></font>
-rw-rw-r-- 1 vidal vidal 1,7K may 13 21:25 README-blastn
-rw-rw-r-- 1 vidal vidal 9,7K may 13 21:25 README.md
-rw-rw-r-- 1 vidal vidal 4,9K may 13 21:25 README-vdb-config
drwxrwxr-x 8 vidal vidal 4,0K may 13 21:32 <font color="#729FCF"><b>schema</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/sratoolkit.3.1.1-ubuntu64</b></font>$ cd ..
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ ls
<font color="#729FCF"><b>anaconda3</b></font>  <font color="#729FCF"><b>Descargas</b></font>  <font color="#729FCF"><b>Documentos</b></font>  <font color="#729FCF"><b>Escritorio</b></font>  <font color="#729FCF"><b>Imágenes</b></font>  <font color="#729FCF"><b>Música</b></font>  <font color="#729FCF"><b>Plantillas</b></font>  <font color="#729FCF"><b>Público</b></font>  <font color="#729FCF"><b>sratoolkit.3.1.1-ubuntu64</b></font>  <font color="#8AE234"><b>stk.tar.gz</b></font>  <font color="#729FCF"><b>Vídeos</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ cd Descargas
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Descargas</b></font>$ ls
<font color="#8AE234"><b>Anaconda3-2024.06-1-Linux-x86_64.sh</b></font>  <font color="#EF2929"><b>google-chrome-stable_current_amd64.deb</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Descargas</b></font>$ cd ..
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~</b></font>$ cd Documentos
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos</b></font>$ ls
<font color="#729FCF"><b>UNMSM</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos</b></font>$ cd unmsm
bash: cd: unmsm: No existe el archivo o el directorio
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos</b></font>$ cd UNMSM
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM</b></font>$ LS
LS: no se encontró la orden
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM</b></font>$ ls
<font color="#729FCF"><b>2024II</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM</b></font>$ cd 2024II
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II</b></font>$ ls
<font color="#729FCF"><b>G_evol</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II</b></font>$ cd G_evol
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol</b></font>$ ls
<font color="#729FCF"><b>Clase2</b></font>
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol</b></font>$ cd Clase2
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol/Clase2</b></font>$ ls
 clase2.txt      commands_2.sh                <font color="#729FCF"><b>ejemplo_1</b></font>   <font color="#729FCF"><b>ERR12389866</b></font>   sra_accessions_1.txt
 commands_1.sh  &apos;Documento sin título.docx&apos;   <font color="#729FCF"><b>ejemplo_2</b></font>   <font color="#729FCF"><b>ERR12543675</b></font>   sra_accessions_2.txt
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol/Clase2</b></font>$ more sra_accessions_1
more: no se puede abrir sra_accessions_1: No existe el archivo o el directorio
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol/Clase2</b></font>$ more sra_accessions_1.txt
ERR12389866
ERR12543675
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol/Clase2</b></font>$ cd ERR12389866/
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol/Clase2/ERR12389866</b></font>$ LS
LS: no se encontró la orden
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol/Clase2/ERR12389866</b></font>$ ls
ERR12389866.sra
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol/Clase2/ERR12389866</b></font>$ mv */*.sra .cd ..
mv: no se puede efectuar `stat&apos; sobre &apos;*/*.sra&apos;: No existe el archivo o el directorio
mv: no se puede efectuar `stat&apos; sobre &apos;.cd&apos;: No existe el archivo o el directorio
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol/Clase2/ERR12389866</b></font>$ mv */*.sra .cd ..
mv: no se puede efectuar `stat&apos; sobre &apos;*/*.sra&apos;: No existe el archivo o el directorio
mv: no se puede efectuar `stat&apos; sobre &apos;.cd&apos;: No existe el archivo o el directorio
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol/Clase2/ERR12389866</b></font>$ cd ..
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol/Clase2</b></font>$ mv */*.sra .cd ..
mv: no se puede efectuar `stat&apos; sobre &apos;.cd&apos;: No existe el archivo o el directorio
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol/Clase2</b></font>$ cd ..
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol</b></font>$ ls
<font color="#729FCF"><b>Clase2</b></font>  ERR12389866.sra  ERR12543675.sra
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol</b></font>$ fasterq-dump --split-files *.sra ;
spots read      : 31,418
reads read      : 61,865
reads written   : 61,865
spots read      : 204,783
reads read      : 409,566
reads written   : 409,566
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol</b></font>$ ls
<font color="#729FCF"><b>Clase2</b></font>  ERR12389866_1.fastq  ERR12389866_2.fastq  ERR12389866.fastq  ERR12389866.sra  ERR12543675_1.fastq  ERR12543675_2.fastq  ERR12543675.sra
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol</b></font>$ ls -lh
total 225M
drwxr-xr-x 6 vidal vidal 4,0K set 15 19:06 <font color="#729FCF"><b>Clase2</b></font>
-rw-rw-r-- 1 vidal vidal 9,2M set 15 19:16 ERR12389866_1.fastq
-rw-rw-r-- 1 vidal vidal 9,2M set 15 19:16 ERR12389866_2.fastq
-rw-rw-r-- 1 vidal vidal 301K set 15 19:16 ERR12389866.fastq
-rw-rw-r-- 1 vidal vidal 3,1M set 15 19:06 ERR12389866.sra
-rw-rw-r-- 1 vidal vidal  88M set 15 19:16 ERR12543675_1.fastq
-rw-rw-r-- 1 vidal vidal  88M set 15 19:16 ERR12543675_2.fastq
-rw-rw-r-- 1 vidal vidal  29M set 15 19:06 ERR12543675.sra
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol</b></font>$ gzip *fastq
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol</b></font>$ ls -lh
total 71M
drwxr-xr-x 6 vidal vidal 4,0K set 15 19:06 <font color="#729FCF"><b>Clase2</b></font>
-rw-rw-r-- 1 vidal vidal 1,6M set 15 19:16 <font color="#EF2929"><b>ERR12389866_1.fastq.gz</b></font>
-rw-rw-r-- 1 vidal vidal 2,6M set 15 19:16 <font color="#EF2929"><b>ERR12389866_2.fastq.gz</b></font>
-rw-rw-r-- 1 vidal vidal  59K set 15 19:16 <font color="#EF2929"><b>ERR12389866.fastq.gz</b></font>
-rw-rw-r-- 1 vidal vidal 3,1M set 15 19:06 ERR12389866.sra
-rw-rw-r-- 1 vidal vidal  17M set 15 19:16 <font color="#EF2929"><b>ERR12543675_1.fastq.gz</b></font>
-rw-rw-r-- 1 vidal vidal  20M set 15 19:16 <font color="#EF2929"><b>ERR12543675_2.fastq.gz</b></font>
-rw-rw-r-- 1 vidal vidal  29M set 15 19:06 ERR12543675.sra
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol</b></font>$ fastqc *
null
Failed to process Clase2
java.io.FileNotFoundException: Clase2 (Es un directorio)
	at java.base/java.io.FileInputStream.open0(Native Method)
	at java.base/java.io.FileInputStream.open(FileInputStream.java:219)
	at java.base/java.io.FileInputStream.&lt;init&gt;(FileInputStream.java:157)
	at uk.ac.babraham.FastQC.Sequence.FastQFile.&lt;init&gt;(FastQFile.java:77)
	at uk.ac.babraham.FastQC.Sequence.SequenceFactory.getSequenceFile(SequenceFactory.java:106)
	at uk.ac.babraham.FastQC.Sequence.SequenceFactory.getSequenceFile(SequenceFactory.java:62)
	at uk.ac.babraham.FastQC.Analysis.OfflineRunner.processFile(OfflineRunner.java:163)
	at uk.ac.babraham.FastQC.Analysis.OfflineRunner.&lt;init&gt;(OfflineRunner.java:125)
	at uk.ac.babraham.FastQC.FastQCApplication.main(FastQCApplication.java:316)
application/gzip
application/gzip
Started analysis of ERR12389866_1.fastq.gz
application/gzip
null
Failed to process ERR12389866.sra
uk.ac.babraham.FastQC.Sequence.SequenceFormatException: ID line didn&apos;t start with &apos;@&apos; at line 1
	at uk.ac.babraham.FastQC.Sequence.FastQFile.readNext(FastQFile.java:163)
	at uk.ac.babraham.FastQC.Sequence.FastQFile.&lt;init&gt;(FastQFile.java:93)
	at uk.ac.babraham.FastQC.Sequence.SequenceFactory.getSequenceFile(SequenceFactory.java:106)
	at uk.ac.babraham.FastQC.Sequence.SequenceFactory.getSequenceFile(SequenceFactory.java:62)
	at uk.ac.babraham.FastQC.Analysis.OfflineRunner.processFile(OfflineRunner.java:163)
	at uk.ac.babraham.FastQC.Analysis.OfflineRunner.&lt;init&gt;(OfflineRunner.java:125)
	at uk.ac.babraham.FastQC.FastQCApplication.main(FastQCApplication.java:316)
application/gzip
application/gzip
null
Failed to process ERR12543675.sra
uk.ac.babraham.FastQC.Sequence.SequenceFormatException: ID line didn&apos;t start with &apos;@&apos; at line 1
	at uk.ac.babraham.FastQC.Sequence.FastQFile.readNext(FastQFile.java:163)
	at uk.ac.babraham.FastQC.Sequence.FastQFile.&lt;init&gt;(FastQFile.java:93)
	at uk.ac.babraham.FastQC.Sequence.SequenceFactory.getSequenceFile(SequenceFactory.java:106)
	at uk.ac.babraham.FastQC.Sequence.SequenceFactory.getSequenceFile(SequenceFactory.java:62)
	at uk.ac.babraham.FastQC.Analysis.OfflineRunner.processFile(OfflineRunner.java:163)
	at uk.ac.babraham.FastQC.Analysis.OfflineRunner.&lt;init&gt;(OfflineRunner.java:125)
	at uk.ac.babraham.FastQC.FastQCApplication.main(FastQCApplication.java:316)
Approx 5% complete for ERR12389866_1.fastq.gz
Approx 10% complete for ERR12389866_1.fastq.gz
Approx 15% complete for ERR12389866_1.fastq.gz
Approx 20% complete for ERR12389866_1.fastq.gz
Approx 25% complete for ERR12389866_1.fastq.gz
Approx 30% complete for ERR12389866_1.fastq.gz
Approx 35% complete for ERR12389866_1.fastq.gz
Approx 40% complete for ERR12389866_1.fastq.gz
Approx 45% complete for ERR12389866_1.fastq.gz
Approx 50% complete for ERR12389866_1.fastq.gz
Approx 55% complete for ERR12389866_1.fastq.gz
Approx 60% complete for ERR12389866_1.fastq.gz
Approx 65% complete for ERR12389866_1.fastq.gz
Approx 70% complete for ERR12389866_1.fastq.gz
Approx 75% complete for ERR12389866_1.fastq.gz
Approx 80% complete for ERR12389866_1.fastq.gz
Approx 85% complete for ERR12389866_1.fastq.gz
Approx 90% complete for ERR12389866_1.fastq.gz
Approx 95% complete for ERR12389866_1.fastq.gz
Analysis complete for ERR12389866_1.fastq.gz
Started analysis of ERR12389866_2.fastq.gz
Approx 5% complete for ERR12389866_2.fastq.gz
Approx 10% complete for ERR12389866_2.fastq.gz
Approx 15% complete for ERR12389866_2.fastq.gz
Approx 20% complete for ERR12389866_2.fastq.gz
Approx 25% complete for ERR12389866_2.fastq.gz
Approx 30% complete for ERR12389866_2.fastq.gz
Approx 35% complete for ERR12389866_2.fastq.gz
Approx 40% complete for ERR12389866_2.fastq.gz
Approx 45% complete for ERR12389866_2.fastq.gz
Approx 50% complete for ERR12389866_2.fastq.gz
Approx 55% complete for ERR12389866_2.fastq.gz
Approx 60% complete for ERR12389866_2.fastq.gz
Approx 65% complete for ERR12389866_2.fastq.gz
Approx 70% complete for ERR12389866_2.fastq.gz
Approx 75% complete for ERR12389866_2.fastq.gz
Approx 80% complete for ERR12389866_2.fastq.gz
Approx 85% complete for ERR12389866_2.fastq.gz
Approx 90% complete for ERR12389866_2.fastq.gz
Approx 95% complete for ERR12389866_2.fastq.gz
Analysis complete for ERR12389866_2.fastq.gz
Started analysis of ERR12389866.fastq.gz
Analysis complete for ERR12389866.fastq.gz
Started analysis of ERR12543675_1.fastq.gz
Approx 5% complete for ERR12543675_1.fastq.gz
Approx 10% complete for ERR12543675_1.fastq.gz
Approx 15% complete for ERR12543675_1.fastq.gz
Approx 20% complete for ERR12543675_1.fastq.gz
Approx 25% complete for ERR12543675_1.fastq.gz
Approx 30% complete for ERR12543675_1.fastq.gz
Approx 35% complete for ERR12543675_1.fastq.gz
Approx 40% complete for ERR12543675_1.fastq.gz
Approx 45% complete for ERR12543675_1.fastq.gz
Approx 50% complete for ERR12543675_1.fastq.gz
Approx 55% complete for ERR12543675_1.fastq.gz
Approx 60% complete for ERR12543675_1.fastq.gz
Approx 65% complete for ERR12543675_1.fastq.gz
Approx 70% complete for ERR12543675_1.fastq.gz
Approx 75% complete for ERR12543675_1.fastq.gz
Approx 80% complete for ERR12543675_1.fastq.gz
Approx 85% complete for ERR12543675_1.fastq.gz
Approx 90% complete for ERR12543675_1.fastq.gz
Approx 95% complete for ERR12543675_1.fastq.gz
Analysis complete for ERR12543675_1.fastq.gz
Started analysis of ERR12543675_2.fastq.gz
Approx 5% complete for ERR12543675_2.fastq.gz
Approx 10% complete for ERR12543675_2.fastq.gz
Approx 15% complete for ERR12543675_2.fastq.gz
Approx 20% complete for ERR12543675_2.fastq.gz
Approx 25% complete for ERR12543675_2.fastq.gz
Approx 30% complete for ERR12543675_2.fastq.gz
Approx 35% complete for ERR12543675_2.fastq.gz
Approx 40% complete for ERR12543675_2.fastq.gz
Approx 45% complete for ERR12543675_2.fastq.gz
Approx 50% complete for ERR12543675_2.fastq.gz
Approx 55% complete for ERR12543675_2.fastq.gz
Approx 60% complete for ERR12543675_2.fastq.gz
Approx 65% complete for ERR12543675_2.fastq.gz
Approx 70% complete for ERR12543675_2.fastq.gz
Approx 75% complete for ERR12543675_2.fastq.gz
Approx 80% complete for ERR12543675_2.fastq.gz
Approx 85% complete for ERR12543675_2.fastq.gz
Approx 90% complete for ERR12543675_2.fastq.gz
Approx 95% complete for ERR12543675_2.fastq.gz
Analysis complete for ERR12543675_2.fastq.gz
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol</b></font>$ ls -lh
total 76M
drwxr-xr-x 6 vidal vidal 4,0K set 15 19:06 <font color="#729FCF"><b>Clase2</b></font>
-rw-rw-r-- 1 vidal vidal 638K set 15 19:18 ERR12389866_1_fastqc.html
-rw-rw-r-- 1 vidal vidal 385K set 15 19:18 <font color="#EF2929"><b>ERR12389866_1_fastqc.zip</b></font>
-rw-rw-r-- 1 vidal vidal 1,6M set 15 19:16 <font color="#EF2929"><b>ERR12389866_1.fastq.gz</b></font>
-rw-rw-r-- 1 vidal vidal 629K set 15 19:18 ERR12389866_2_fastqc.html
-rw-rw-r-- 1 vidal vidal 407K set 15 19:18 <font color="#EF2929"><b>ERR12389866_2_fastqc.zip</b></font>
-rw-rw-r-- 1 vidal vidal 2,6M set 15 19:16 <font color="#EF2929"><b>ERR12389866_2.fastq.gz</b></font>
-rw-rw-r-- 1 vidal vidal 768K set 15 19:18 ERR12389866_fastqc.html
-rw-rw-r-- 1 vidal vidal 458K set 15 19:18 <font color="#EF2929"><b>ERR12389866_fastqc.zip</b></font>
-rw-rw-r-- 1 vidal vidal  59K set 15 19:16 <font color="#EF2929"><b>ERR12389866.fastq.gz</b></font>
-rw-rw-r-- 1 vidal vidal 3,1M set 15 19:06 ERR12389866.sra
-rw-rw-r-- 1 vidal vidal 578K set 15 19:18 ERR12543675_1_fastqc.html
-rw-rw-r-- 1 vidal vidal 325K set 15 19:18 <font color="#EF2929"><b>ERR12543675_1_fastqc.zip</b></font>
-rw-rw-r-- 1 vidal vidal  17M set 15 19:16 <font color="#EF2929"><b>ERR12543675_1.fastq.gz</b></font>
-rw-rw-r-- 1 vidal vidal 579K set 15 19:18 ERR12543675_2_fastqc.html
-rw-rw-r-- 1 vidal vidal 329K set 15 19:18 <font color="#EF2929"><b>ERR12543675_2_fastqc.zip</b></font>
-rw-rw-r-- 1 vidal vidal  20M set 15 19:16 <font color="#EF2929"><b>ERR12543675_2.fastq.gz</b></font>
-rw-rw-r-- 1 vidal vidal  29M set 15 19:06 ERR12543675.sra
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol</b></font>$ fastqc
application/gzip
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol</b></font>$ fastqc
application/gzip
(base) <font color="#8AE234"><b>vidal@PC</b></font>:<font color="#729FCF"><b>~/Documentos/UNMSM/2024II/G_evol</b></font>$ 
</pre>
