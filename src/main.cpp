// arpackcpp.cpp : Defines the entry point for the console application.
//


#include "stdafx.h"
//#include "D:\ap\arpack++\examples\matprod\nonsym\nmatrixv.h"
#include "D:\ap\arpack++\examples\reverse\nonsym\rnsymvsl.h"
#include "arrsnsym.h"
#include <vector>
#include <ctime>
#include <string>
#include<iostream>

using namespace std;


int numRows;
int numColumns;
int numNonZeroValues;

char *precision = "SINGLE";
char *whichVectors = "LM";
char *type = "SVDLEFT";


class eigValVectorPair
{   
   public: 
      double magnitude;
	  int vectorColumn;
};


// multiplies matrix M^T with v
template<typename T>
inline void MultMT(T* v, T* w, vector<long>& rowsptr, vector<long>& pointr, vector<T>& values)
{ 
  int i,j;
  for(i=0;i<numColumns;i++) w[i]=0;
  //pomnozimo sa A'  (nxm) puta (mx1) ispadne (nx1)
  for(i=0;i<numColumns;i++){ //za svaki stupac
	  for(j=pointr[i];j<pointr[i+1];j++) { //sumiramo elemente u stupcu pomnozene odgovarajucim elementima od v
		w[i]+=values[j]*v[rowsptr[j]];
	  }
  }
}

// multiplies matrix M with v
template<typename T>
inline void MultM(T* v, T* w, vector<long>& rowsptr, vector<long>& pointr, vector<T>& values)

{ 
  int i,j;
  for(i=0;i<numRows;i++) w[i]=0;
  for(i=0;i<numColumns;i++){ //za svaki stupac
	  for(j=pointr[i];j<pointr[i+1];j++){
		w[rowsptr[j]]+=values[j]*v[i];
	  }
  }
} // MultOPv


template<typename T>
inline void MultOPv(T* v, T* w, vector<long>& rowsptr, vector<long>& pointr, vector<T>& values) 
/*
  Performs the matrix-vector product w <- X*v, with X being matrix M or MM' or M'M (M' is the transpose of M)
*/
{   
  if(strcmp(type,"EIGDEC")==0){ // eigenvectors of M
	  MultM(v,w,rowsptr,pointr,values);
  } else if(strcmp(type,"SVDLEFT")==0){ // eigenvectors of MM'
	T* t = new T[numColumns]();
	MultMT(v,t,rowsptr,pointr,values);
	MultM(t,w,rowsptr,pointr,values);
	delete[] t;
  } else if(strcmp(type,"SVDRIGHT")==0){ // eigenvectors of M'M
	T* t = new T[numRows]();
	MultM(v,t,rowsptr,pointr,values);
	MultMT(t,w,rowsptr,pointr,values);
	delete []t;
  }
} // MultOPv


template<typename T> //loads data from the file, runs the arpack routines and writes the results
void Solve(char* inFileName, char* outFileName, long numDim)
{

	
	FILE *infile = fopen(inFileName,"rb");
	FILE *outf = fopen(outFileName,"w");

	if(!infile){
		printf("Can't open file: %s", inFileName);
		return;
	} 

	if(!outf){
		printf("Can't open file: %s",outFileName);
		return;
	}

	cout<<"Loading dimensions."<<endl;
	
	fread(&numRows,sizeof(int),1,infile);
	fread(&numColumns,sizeof(int),1,infile);
	fread(&numNonZeroValues,sizeof(int),1,infile);

	cout<<"Initialising rowsptr..."<<endl;
	vector<long> rowsptr(numNonZeroValues);
	cout<<"Initialising pointr..."<<endl;
	vector<long>pointr(numColumns + 1);
	cout<<"Initialising values ..."<<endl;
	vector<T> values(numNonZeroValues);
  
	cout<<"Loading values ..."<<endl;	
	int intvar;
    T valuevar;
	for(int i = 0; i < numNonZeroValues; i++)
	{
		//fgets(str,50,infile);
		if(i%100000 == 0) cout<<"Value "<<i<<" of "<<numNonZeroValues<<endl;		
		//fscanf(infile,"%s",str);
		//values[i] = atof(str);
		
		fread(&valuevar,sizeof(T),1,infile);
		values[i]=valuevar;
		//cout<<values[i]<<endl;
	}
	cout<<" done."<<endl;

	cout<<"Loading rows ..."<<endl;
	for(int i = 0; i < numNonZeroValues; i++)
	{
		if(i%100000 == 0) cout<<"Row"<<i<<" of "<<numNonZeroValues<<endl;
		//fscanf(infile,"%s",str);
		//rowsptr[i]=atol(str);		
		fread(&intvar,sizeof(int),1,infile);
		rowsptr[i]=intvar;
	}
	cout<<" done."<<endl;

	cout<<"Loading columns ... "<<endl;
	for(int i = 0; i <= numColumns; i++)
	{
		if(i%100000 == 0) cout<<"Column "<<i<<" of "<<numColumns<<endl;
		//fscanf(infile,"%s",str);
		//pointr[i]=atol(str);
		fread(&intvar,sizeof(int),1,infile);
		pointr[i]=intvar;
	}
	cout<<" done."<<endl;
	        
	fclose(infile);

	// for SVDRIGHT the underlying quadratic matrix for which we need eigenvectors is numColums x numColumns
	// for SVDRIGHT it is numRows x numRows
	// for EIGENDEC numRows must be equal to numColumns 	
	int matDim = strcmp(type,"SVDRIGHT")==0 ? numColumns : numRows; 	
	//some more error checking
	if(strcmp(type,"EIGENDEC")==0 && numColumns!=numRows) throw "EIGENDEC option allows only square matrices.";
    if(!(numDim>1 && numDim<(matDim-1))) throw "numDim must be >1 and <matrixDimension";
	
	//initialize the arpack++ class which provides a nice interface to the underlying fortran routines	    
	
	ARrcNonSymStdEig<T> prob(matDim, numDim);
    prob.ChangeWhich(whichVectors);

    // Finding an Arnoldi basis.
    int iter = 1;
    printf("Starting iterations ...\n");
    while (!prob.ArnoldiBasisFound()) {
    // Calling ARPACK FORTRAN code. Almost all work needed to
    // find an Arnoldi basis is performed by TakeStep.
		prob.TakeStep();
		if ((prob.GetIdo() == 1)||(prob.GetIdo() == -1)) {
			MultOPv(prob.GetVector(), prob.PutVector(),rowsptr,pointr,values);
    }
    printf("Finished iteration: %d\n",iter);
	iter++;
  }
  cout<<"Finished iterations!"<<endl;
  
  // Finding eigenvalues and eigenvectors.
  prob.FindEigenvectors();
  
  // print some information
  Solution(prob);

  // Sort the vectors descending
  int nconv = prob.ConvergedEigenvalues();
  vector<eigValVectorPair> eigenValues(nconv); 
  if (prob.EigenvaluesFound()) {
 
    for (int i=0; i<nconv; i++) {
		eigenValues[i].magnitude = sqrt(prob.EigenvalueReal(i) * prob.EigenvalueReal(i) + prob.EigenvalueImag(i)*prob.EigenvalueImag(i));
		eigenValues[i].vectorColumn = i;
    } 
	for(int i = 0; i < nconv - 1; i++){
		for(int j = i + 1; j < nconv; j++){
			if(eigenValues[i].magnitude < eigenValues[j].magnitude){
				eigValVectorPair tmp = eigenValues[i];
				eigenValues[i] = eigenValues[j];
				eigenValues[j] = tmp;
			}
		}
	}
  }

  // Printing eigenvalues and eigenvectors to the outputfile.  
  printf("Writing vectors to file.\n");  
  fprintf(outf,"%d %d\n",numRows, numDim);
  for(int i = 0; i <numRows; i++){
	  for(int j = 0; j < numDim; j++)
		  fprintf(outf,"%f ",prob.EigenvectorReal(eigenValues[j].vectorColumn,i));
	  fprintf(outf,"\n");
  }
  fclose(outf);

} // Test.


int main(int argc, char* argv[])
{
  clock_t start, finish;
  start = clock();

  char *inFileName, *outFileName;
  long d;
  //parsing input arguments 
  try{	  
	 inFileName = argv[1];
	  outFileName  = argv[2];
	  d = atol(argv[3]);
	  printf("d je %d", d);
	  for(int i = 4;i<argc;i+=2){
		  if(strcmp(argv[i], "-t")==0){
			if(strcmp(argv[i+1],"EIGDEC")==0 || strcmp(argv[i+1],"SVDLEFT")==0 || strcmp(argv[i+1],"SVDRIGHT")==0)
				type = argv[i+1];
			else throw "Type of decomposition not supported (only allowed are EIGDEC, SVDLEFT and SVDRIGHT";
		  } else if(strcmp(argv[i],"-p")==0){
			if(strcmp(argv[i+1],"SINGLE")==0 || strcmp(argv[i+1],"DOUBLE")==0)
				precision = argv[i+1];
			else throw "Precision not supported (only allowed options are SINGLE and DOUBLE";
		  } else if(strcmp(argv[i],"-w")==0){
			if(strcmp(argv[i+1],"LM")==0 || strcmp(argv[i+1],"SM")==0)
				whichVectors = argv[i+1];
			else throw "WHICH option not supported (only allowed options are LM and SM)";
		  }
	  }
	  //inFileName = "C:\\xlansimdata\\LSAMatrices\\europarl-v7.de-en.en.ende.lsamatrix.bin";
	  //outFileName = "C:\\xlansimdata\\LSADecomps\\test2.txt";

	  if(strcmp(precision,"SINGLE")==0)
		 Solve<float>(inFileName, outFileName,d);
	  else
		 Solve<double>(inFileName, outFileName,d);
  } catch (const std::exception& ex) {
	  cout << "Exception raised: " << ex.what() << endl;
  } catch (const std::string& ex) {
	  cout << "Exception raised: " << ex << endl;
  } catch (...) {
    cout << "Unknown exception raised!" << endl;
	return -2;
  }


  
  printf("All done ^_^."); 
  finish = clock();
  float sec = ((finish - start)/CLOCKS_PER_SEC );
  printf("\nTime taken: %f seconds - %f minutes",sec,sec/60.0);

} // main

