Spectral-Unmixing
=================

Algorithms for spectral unmixing in C programming language using OpenMP 


|------------------Instructions------------------|



In order to use the programs you must have GDAL installed. GDAL is used to convert the input image files to .ers extention, so that the programs can read them. 


The programs use comand line arguments. Use them by entering the corresponding path to the terminal, so the program can inform you about the arguments it needs. After that enter the arguments following the path to the program.

 


|--------------Compiler Instructions--------------|




GNU/Linux 

 

1) Open a terminal in the folder where you put the code files.

2) Run: gcc code_file_name.c -o executable_file name -lm -fopenmp -Ofast

 

Microsoft Windows 

 

Compiler Instructions are referred to the open source IDE CodeBlocks, available at www.codeblocks.org 

 

1) Download and install TDM-GCC checing the correcponding openMP box.

3) In Codeblocks, choose Settings->Compiler. 

4) Choose GNU GCC compiler and then copy and rename. 

5) At linker settings, at link libraries choose add point to libgomp-1.dll which is in the bin directory of TDM-GCC folder. 

6) At Toolchain executables, at Compiler's installation Directory write the TDM-GCC path. 

7) At compiler settings, at other settings write "-fopenmp" and "-Ofast". 




|--------------------Text Files-------------------|

 

Text files used as input files must follow a specific patern. The first line contains user data, the programs ignore it. In the next, there must be the value of all the endmembers in each band. The values of each endmember in a specific band are space delimited. So, each endmember occupies a column. For Example if we have 5 endmembers : 


Endmembers : 

245.55 154.78 16.25 165.13 12.74

214.68 164.35 23.44 143.87 20.46

.        .     .      .     .

.        .     .      .     .

.        .     .      .     .
.        .     .      .     .

194.74 241.98 11.20 182.41 14.08


|--------------------Programs-------------------|


--ODM.c : Performs the Outlier Detection Method (ODM) in order to estimate the number of endmembers of a hyperspectral image.

--SEE.c : Performs the Simple Endmember Extraction Method (SEE) in order to extract the spectral signatures of the endmembers of a hyperspectral image.

--SEE-E.c : Performs the Simple Endmember Extraction Method -Enchanced (SEE-E) in order to extract the spectral signatures of the endmembers of a hyperspectral image.

--AFEM.c : Performs Abundance Fraction Estimation Method in order to estimate the abundance fraction of each endmember in every pixel of a hyperspectral image.

--MNF.c : Performs the Maximum Noise Fraction Method (MNF) on a hyrespectral image.

--PCA.c : Performs the Principal Components Transformation (PC) on an Image.

--Unmixing.c : Performs MNF, ODM, SEE-E and AFEM Methods on a hyrespectral image.

--SAD.c : Creates an image where the value every pixel is the spectral angle of its true value with a given set of signatures.

--SADtxt.c : Calculates the spectral angles of two sets of signatures.
