The wsvd tool was made to perform efficient computation of truncated SVD as 
efficiently as possible for very large sparse matrices. It is a wrapper around 
the svd routines from the ARPACK library.

Author: Mladen Karan, Text Analysis and Knowledge Engineering Lab, Faculty of Electrical Engineering and Computing, University of Zagreb

1. Usage:
wsvd input output K [-t TYPE] [-w WHICH] [-p PREC]

input - the input file name
output - the output file name
K - number of vectors to compute
t - what type of vectors is desired TYPE can be (default is SVDLEFT):
     - SVDLEFT - left singular vectors
     - SVDRIGHT - right singular vectors
     - EIGDEC - computes eigenvectors
w - which vectors are desired WHICH can be (default is LM)
      - LM - those corresponding to the largest magnitude singular values
      - SM - those corresponding to the smallest magnitude singular values
p - precision at which to operate (default is SINGLE)
      - SINGLE - single precision (float)
      - DOUBLE - double precision (double)





2. Input file formats:

a)

The wsvd utility itself works with the input files in compressed sparse column (CSC) format.
However, to make it easier to work with word/context matrices the following format is also possible.

word1 columnIndex:value columnIndex:value ...
word2 columnIndex:value columnIndex:value ...

Examples of such matrices are in the examples subfolder.
Before using wsvd on such a matrix it must be converted to the binary format. A utility for this is called conv2CSC.exe 
and is provided in the package. A matrix mat.txt with a total of numCol columns can be converted to mat.bin as follows:

conv2CSC.exe mat.txt mat.bin numOfColumns

b)

The binary format itself is as follows (for a matrix with r rows, c columns and nz non-zero elements):

1 int - number of rows (r)
1 int - number of columns (c)
1 int - number of non-zeros (z)
nz float/double - the values (the value array)
nz int - the rows corresponding to the values (the rowindexes array)
c+1 int - for each column the index in the previous two arrays where that column starts (including an entry for a fictional c+1 column)

See these links for more details on the CSC format:
http://netlib.org/linalg/html_templates/node91.html
http://netlib.org/linalg/html_templates/node92.html




3. Output file format
The output file contains the desired vectors as follows:

ndims nvectors
singval1
singval2
singval3
...
v1component1   v2component1    v3component1 ...
v1component2   v2component2    v3component2 ...
...	       ... 	       ...

The vectors are sorted in order of decreasing magnitude of the corresponding singular values.
  



4. Example - LSA
Assume we have some words and contexts in a sparse matrix (in examples\exampleSmall.txt, it has 8 columns)

truck 1:9.0	8:15
panda 7:13
shroud 4:12	8:7
phone 1:3
wall 2:16
key 5:5
violin 6:4
space 3:2

1. we convert it to the binary format -> convCSC example\exampleSmall.txt exampleSmall.bin 8
2. and perform the svd (finding 3 most important vectors) -> wsvd exampleSmall.bin output.txt 3
   - we want the largest magnitude left singular vectors using single precision - these are exactly the default settings


If everything went well the output file should look like this
8 3
19.241273
16.000001
13.000001
-0.858567 -0.000000 -0.000000 
-0.000000 -0.000000 1.000000 
-0.508668 0.000000 0.000000 
-0.064174 0.000000 0.000000 
-0.000000 -1.000000 -0.000000 
-0.000000 0.000000 0.000000 
0.000000 0.000000 0.000000 
0.000000 0.000000 -0.000000 

The 3 left singular vectors make the truncated U matrix from SVD and we can 
read the rows of this matrix as the word representations. So the LSA vector
corresponding to shroud would be [-0.508668 0.000000 0.000000] 



5. Additional remarks
 	- for the binary format indexes are 0 based, for the word/context format they are 1 based
	- conv2CSC will interpret values as floats
	- the non LSA specific settings (e.g. getting the right singular vectors) have not been tested
          well and it is recommended to double check correctness of results on a smaller matrix using 
          other software (e.g. R) before running on full matrices. 

Please report any bugs to mladen.karan@fer.hr



