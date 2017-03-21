// Threaded two-dimensional Discrete FFT transform
// Author:Amit Sunil Dhamne
// ECE8893 Project 2


#include <iostream>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "Complex.h"
#include "InputImage.h"

// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.

using namespace std;

Complex buf[1024*1024];
Complex* h; Complex* h1; 
Complex W[512];
Complex H1[1024*1024];
Complex H2[1024*1024]; 
Complex Wi[512]; 
int r,c;
int p,count;
bool* localSense;
bool globalSense;

int FetchAndDecrementCount();
void Transpose(Complex* h, int N);
int counts=0;
int nThreads = 16;
// The mutex and condition variables allow the main thread to
// know when all helper threads are completed.
pthread_mutex_t transposeMutex;
pthread_mutex_t threadCount;
pthread_mutex_t countMutex;
pthread_mutex_t exitMutex;
pthread_cond_t  exitCond;
 

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.
unsigned ReverseBits(unsigned int v, int N)
{ //  Provided to students
  unsigned n = N; // Size of array (which is even 2 power k value)
  unsigned r1 = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r1 <<= 1;        // Shift return value
      r1 |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r1;
}

// GRAD Students implement the following 2 functions.
// Undergrads can use the built-in barriers in pthreads.

// Call MyBarrier_Init once in main
void MyBarrier_Init()// you will likely need some parameters)
{
  p = nThreads;
  count = nThreads;
  pthread_mutex_init(&countMutex,0);
  localSense = new bool[p];
  for(int i=0;i<p;i++) localSense[i]=true;
  globalSense = true;
}

// Each thread calls MyBarrier after completing the row-wise DFT
void MyBarrier(int myId) // Again likely need parameters
{
  localSense[myId] = !localSense[myId];
  if(FetchAndDecrementCount() == 1)
    {
      count =p;
      globalSense = localSense[myId];
    }
  else
    {
      while (globalSense != localSense[myId]) {}
    }
}
int FetchAndDecrementCount()
{
  pthread_mutex_lock(&countMutex);
  int myCount = count;
  count--;
  pthread_mutex_unlock(&countMutex);
  return myCount;
}
       
void Transform1D(Complex* h, int N)//, int i)//, int i)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)
  Complex tempi,tempr;
  Complex temp; int tl;
  Complex H[N];
  for(int y=0; y<=N;y++) 
    {
      tl=  ReverseBits(y,N);
      H[y] = h[tl];  
    }
   for(int j=2;j<=N;j=j*2)
     {
      for(int jj=0;jj<N/j;jj = jj++)
        {
          for(int k =jj*j ; k<(jj*j)+j/2;k=k++)
             {  
	        tempi = H[k] + W[N*(k-jj*j)/j]*H[k+j/2];
	        tempr  = H[k] - W[N*(k-jj*j)/j]*H[k+j/2];
          	H[k] = tempi; 
          	H[k+j/2]=tempr;
              }
          }
       }
   copy(H,H+N,h);
   
}
void Inverse(Complex*h, int N)
{
  Complex tempi, tempr;
  int tl;Complex temp;
  temp.imag,temp.real=1.0/N;
  Complex H[N];
  for(int x=0;x<N;x++)
    {
      tl= ReverseBits(x,N);
     H[x] = h[tl];
    }
  for(int j=2;j<=N;j=j*2)
    {                                                                                                                                                                          
      for(int jj=0;jj<N/j;jj = jj++)
        {                                                                                                                                                                      
          for(int k =jj*j ; k<(jj*j)+j/2;k=k++)
	    {                                                                                                                                                                   
	      tempi = H[k] + Wi[N*(k-jj*j)/j]*H[k+j/2];
	      tempr  = H[k] - Wi[N*(k-jj*j)/j]*H[k+j/2];
	      H[k] = tempi;
	      H[k+j/2]=tempr;
	    }
	}
    }
  for (int u=0;u<N;u++) H[u]=H[u]*temp;
  copy(H,H+N,h);
    
}

void* Transform2DTHread(void* v)
{ // This is the thread startign point.  "v" is the thread number
  // Calculate 1d DFT for assigned rows
  // wait for all to complete
  // Calculate 1d DFT for assigned columns
  // Decrement active count and signal main if all complete
  
  unsigned long myId =  (unsigned long)v;
  int rowsPerThread = c/nThreads;
  int startingRow = myId*rowsPerThread;
  for(int i=0;i<rowsPerThread;i++)
    {
      Transform1D(&h[(startingRow+i)*c],c);
    } 
  MyBarrier(myId);
  if (myId==0)
    {
      copy(h,h+c*c,H1);h=h1;
      pthread_mutex_lock(&transposeMutex);
      Transpose(h,c);
      pthread_mutex_unlock(&transposeMutex);
    }
  MyBarrier(myId);
  for(int j=0;j<rowsPerThread;j++)
     {
       Transform1D(&h[(startingRow+j)*c],c);
     }
  MyBarrier(myId);
  if(myId==0)
    {
      pthread_mutex_lock(&transposeMutex);
      Transpose(h,c);
      pthread_mutex_unlock(&transposeMutex);
      copy(h,h+c*c,H2);h=h1;
    }
  MyBarrier(myId);
  for(int x = 0;x<rowsPerThread;x++)
    {
      Inverse(&h[(startingRow+x)*c],c);
    }
  MyBarrier(myId);
  if(myId==0)
    {
      pthread_mutex_lock(&transposeMutex);
      Transpose(h,c);
      pthread_mutex_unlock(&transposeMutex);
    } MyBarrier(myId);
  for(int y =0;y<rowsPerThread;y++ )
    {
      Inverse(&h[(startingRow+y)*c],c);
    } MyBarrier(myId);
  if(myId==0)
    {
      pthread_mutex_lock(&transposeMutex);
      Transpose(h,c);
      pthread_mutex_unlock(&transposeMutex);
    } MyBarrier(myId);
  
  pthread_mutex_lock(&threadCount);
  counts++;
  if(counts==nThreads)
    { cout<<"Count"<<counts<<endl;
      pthread_mutex_unlock(&threadCount);
      pthread_mutex_lock(&exitMutex);
      pthread_cond_signal(&exitCond);
      pthread_mutex_unlock(&exitMutex);
    }
  else
    {pthread_mutex_unlock(&threadCount);}
  return 0;
}


void Transform2D(const char* inputFN) 
{ // Do the 2D transform here.
   InputImage image(inputFN);  // Create the helper object for reading the image
  // Create the global pointer to the image array data
  // Create 16 threads
  // Wait for all threads complete
  // Write the transformed data
   c = image.GetWidth();
   r = image.GetHeight();
   h  = image.GetImageData();
   h1= h;
  cout<<r<<" rows\t"<<c<<" column"<<endl;
  pthread_mutex_init(&exitMutex,0); 
  pthread_cond_init(&exitCond,0);
  pthread_mutex_init(&transposeMutex,0);
  pthread_mutex_init(&threadCount,0);
  pthread_mutex_lock(&exitMutex); 
  for(int j = 0;j<c/2;j++)
     {
       W[j] = Complex(cos(2*M_PI*j/c),-sin(2*M_PI*j/c));
    }
 
  for(int jj = 0;jj<c/2;jj++)
    {
      Wi[jj] = Complex(cos(2*M_PI*jj/c),sin(2*M_PI*jj/c));
    }
  pthread_t pt[nThreads];
  MyBarrier_Init();
  for(int i =0;i<nThreads;i++)
    {
      pthread_create(&pt[i],0,Transform2DTHread,(void*)i);
    }
  pthread_cond_wait(&exitCond,&exitMutex);
                                                                                                                             
     h=h1;     
     image.SaveImageData("MyAfter1D.txt",H1,r,c);
     image.SaveImageData("MyAfter2D.txt",H2,c,r);
     image.SaveImageData("MyAfterInverse.txt",h,c,r);  
}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  // MPI initialization here
  Transform2D(fn.c_str()); // Perform the transform.
}  
void Transpose(Complex* h,int N)
{
      for(int i=0;i<N;i++)
	{
	  for(int j=0; j<N; j++)
	    {
	      buf[j+i*N] = h[j*N+i];
	    }
    
          }
      copy(buf,buf+(N*N),h);
}
  
