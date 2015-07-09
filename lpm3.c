/* Copyright (c) 2015  Jonathan Lisic 
 * Last edit: 15/07/08 - 11:45:28
 * License: BSD 3-Clause
 */  


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#ifdef CLI
#include <time.h>
# else 
#include "R.h"
#include "Rmath.h"
#endif

#include "kdtree_lpm.h"


/* return new prob */
void updateProb( double * xPtr, double * yPtr, size_t * xIndex, size_t * yIndex, size_t n ) { 

  double xySum = *xPtr + *yPtr;

  #ifdef CLI
  double U = (double) rand()/ (double) RAND_MAX ; 
  #else
  double U;
  U = runif(0.0, 1.0);
  #endif

  if( xySum < 1 ) {
    if( *yPtr / xySum < U ) {
      *xPtr = 0;
      *yPtr = xySum;
      *xIndex = n;
      return;
    } 

    *xPtr = xySum;
    *yPtr = 0;
    *yIndex = n;
    return;
  } 

  if( (1 - *yPtr) / (2 - xySum) < U ) {
    *xPtr = 1;
    *yPtr = xySum - 1;
    *xIndex = n;
    if( xySum == 1) *yIndex = n;
    return;
  } 
    
  *yPtr = 1;
  *xPtr = xySum - 1;
  *yIndex = n;
  if( xySum == 1) *xIndex = n;
  return;
}





#ifdef CLI

/* a simple C test case */
int main ( int argc, char * argv[]  )  {

  size_t n = 100000;
  size_t m = 40; // leaf size
  size_t K = 2; // dimension
  double sampleRate  = 0.10;
  size_t i,k;

  int nInt;
  int KInt;
  int mInt;

  printf("argc = %d\n", argc);

  if( argc == 5) { 
    sscanf(argv[1], "%d", &nInt); 
    sscanf(argv[2], "%d", &KInt); 
    sscanf(argv[3], "%d", &mInt); 
    sscanf(argv[4], "%lf", &sampleRate);  // c99 requirement on lf
    
    n          = (size_t) nInt;
    m          = (size_t) mInt;
    K          = (size_t) KInt; 
    printf(" Parameters n = %d, m = %d, K = %d, pi = %f\n", (int) n, (int) K, (int) m, sampleRate); 
  }


  /* allocate some memory and make a junk data set */ 
  double * x = calloc( n * K, sizeof(double));
  double * pi = calloc( n , sizeof(double));
  
  //srand(0);
  srand( time(NULL) );
  
  /***************************** GENERATE FAKE DATA ****************************/
  for( i = 0; i < n; i ++) {
    //for( k = 0; k < K; k ++) x[i*K + k]    = (double) i;
    for( k = 0; k < K; k ++) x[i*K + k]    = (double) 1;
    pi[i]     = sampleRate;
  }
  /*************************END GENERATE FAKE DATA ****************************/


#else

/* the R interface */
  void R_lpm3(
      double * x,
      double * pi,
      int * nPtr,
      int * KPtr, 
      int * mPtr
      ) {

  size_t n = (size_t) * nPtr;
  size_t m = (size_t) * mPtr;
  size_t K = (size_t) * KPtr;
  size_t i,k;

  GetRNGstate();

#endif



  size_t j,l; // incrementor
  size_t sampled;
  double dist;
  
  size_t numberIndexToRemove;    // number of indexes to remove
  size_t removeIndexA;           // first index to remove
  size_t removeIndexB;           // second index to remove
  size_t removeReduceIndexA;     // first reduced index to change 
  size_t removeReduceIndexB;     // second reduced index to change 
  size_t numberChosenIndexAtEnd; // number of indexes to move
  size_t moveA;                  // first index to move 
  size_t moveB;                  // second index to move



  /***************************** CREATE INDEX ****************************/

  size_t * reduceIndexMap = (size_t *) calloc( n , sizeof( size_t ) );
  size_t * indexMap = (size_t *) calloc( n , sizeof( size_t ) );
  for( i=0; i< n; i++) reduceIndexMap[i]=i, indexMap[i]=i;
  
  /***************************** CREATE INDEX ****************************/



  rootNodePtr myTree = createTree( K, m, n, x);

  buildIndex( myTree, NULL, NULL); 

#ifdef DDEBUG
printf(" --------- Print Tree -----------\n" );
printTree( myTree, myTree->root ); 
#endif 

  /* algorithm */
  for( i = 0; i < n; i++) {

    /* randomly select j */
#ifdef CLI
    sampled = (size_t) (n - i) * ( (double) rand() / (double) RAND_MAX ) ;
#else 
    sampled = (size_t) floor( runif(0.0, (double) (n - i) ) );
#endif
    j = reduceIndexMap[sampled]; // j is the original index of the sampled data

    if( j >=n  ) {
      //printf("breaking on iteration %d, for k\n", (int) i);
      break;
    }

  
#ifdef DDEBUG

printf("----------------------- %d  (%d, %d)---------------\n", (int) i, (int) j, (int) sampled); 

printf(" %d:\tpoint tree\ttree\tindex\treduce\tpi\n", (int) i); 
for( l = 0; l < n; l ++) 
printf(" %d:\t%p\t%d\t%d\t%d\t%f\n", 
(int) l, 
(void *) myTree->pointerIndex[l],
(int) *( myTree->pointerIndex[l] ),
(int) ( indexMap[l] ),
(int) ( reduceIndexMap[l] ),
pi[l]
);
#endif

    // find neighbor
    dist = INFINITY;
    k = n;
    find_nn_notMe( myTree , myTree->root, j, &dist, &k, &x[j*K]); 

#ifdef DDEBUG
printf(" neighbor: %d\n", (int) k );
#endif

    /* break if an invalid neighbor is selected */
    if( k >=n  ) {
      printf("breaking on iteration %d, for k\n", (int) i);
    break;
    }

  


    updateProb( 
       &( pi[k] ), 
       &( pi[j] ), 
       myTree->pointerIndex[k], 
       myTree->pointerIndex[j], 
       n); 

    /*
    
      the point of index is to provide a mapping from the original n items to the current
      position in reduce index.
      Likewise the point of reduce index is to provid a mapping from the reduce set use
      for sampling to the position in the index set, which is also the mapping to x.
    
    
      1. Figure out how many values must be set to n.
    */


    numberIndexToRemove = 1;
    removeIndexA = k;
    removeIndexB = n;
    removeReduceIndexA = indexMap[k]; 
    removeReduceIndexB = n;
  
    if ( (pi[j] <= 0) | (pi[j] >= 1) ) {
      removeIndexA = j;
#ifdef DDEBUG
printf(" sampled: %d\n", (int) sampled);
#endif
      removeReduceIndexA = sampled; 
      if ( (pi[k] <= 0) | (pi[k] >= 1) ) {
        removeIndexB = k;
        removeReduceIndexB = indexMap[k]; 
        numberIndexToRemove = 2;
      }
    }
 
#ifdef DDEBUG
printf("  removeIndexA = %d\t removeReduceIndexA = %d\n", (int) removeIndexA, (int) removeReduceIndexA);
printf("  removeIndexB = %d\t removeReduceIndexB = %d\n", (int) removeIndexB, (int) removeReduceIndexB);
printf("  numberIndexToRemove = %d\n", (int) numberIndexToRemove);
#endif

    /*       
    2. Figure out which values must be moved. 
    */
  
    moveA = n;
    moveB = n;

    if( numberIndexToRemove == 2 ) {
      moveA = reduceIndexMap[n - i - 1];
      moveB = reduceIndexMap[n - i - 2];
    } else {
      moveA = reduceIndexMap[n - i - 1];
    } 

#ifdef DDEBUG
printf("  moveA = %d\t moveB = %d\n", (int) moveA, (int) moveB );
#endif

    /*
      3. If either of these move values are k or j
         we are just going to set them to n.
    */

    numberChosenIndexAtEnd = 0;

    if( numberIndexToRemove == 2 ) {

      if( removeIndexA == moveA ) 
      {
#ifdef DDEBUG
printf("3a\n");
#endif
        reduceIndexMap[ indexMap[moveA] ] = n;
        indexMap[moveA] = n;
        moveA = n;
        numberChosenIndexAtEnd++;
      }
      else if( removeIndexB == moveA ) 
      {
#ifdef DDEBUG
printf("3b\n");
#endif
        reduceIndexMap[ indexMap[moveA] ] = n;
        indexMap[moveA] = n;
        moveA = n;
        numberChosenIndexAtEnd++;
      }
      
      if( removeIndexA == moveB ) 
      {
#ifdef DDEBUG
printf("3c\n");
#endif
        reduceIndexMap[ indexMap[moveB] ] = n;
        indexMap[moveB] = n;
        moveB = n;
        numberChosenIndexAtEnd++;
      } 
      else if( removeIndexB == moveB ) 
      {
#ifdef DDEBUG
printf("3d\n");
#endif
        reduceIndexMap[ indexMap[moveB] ] = n;
        indexMap[moveB] = n;
        moveB = n;
        numberChosenIndexAtEnd++;
      }
      
      
    } else { 

      if( removeIndexA == moveA ) {
#ifdef DDEBUG
printf("3e\n");
#endif
        reduceIndexMap[ indexMap[moveA] ] = n;
        indexMap[moveA] = n;
        moveA = n;
        numberChosenIndexAtEnd++;
      } 
      else if( removeIndexB == moveA ) 
      {
#ifdef DDEBUG
printf("3f\n");
#endif
        reduceIndexMap[ indexMap[moveA] ] = n;
        indexMap[moveA] = n;
        moveA = n;
        numberChosenIndexAtEnd++;
      }

    }

    /*
      4.  Handle moving is the last set, which we only do
          if j or k is less than n - i - numberIndexToRemove
    */

    if( numberIndexToRemove > numberChosenIndexAtEnd ) {
  
      // check if moveA has been set to n
      if( moveA < n ) {
        // check if removeIndexA has been set to n
        if( indexMap[removeIndexA] != n ) {
#ifdef DDEBUG
printf("4a\n");
#endif
          // set reduced index
          reduceIndexMap[ removeReduceIndexA ] = reduceIndexMap[ indexMap[ moveA ] ];
          // set index
          indexMap[ moveA ] = removeReduceIndexA;
          indexMap[ removeIndexA ] =n;
        } else {
#ifdef DDEBUG
printf("4b\n");
#endif
          reduceIndexMap[ removeReduceIndexB ] = reduceIndexMap[ indexMap[ moveA ] ];
          indexMap[ moveA ] = removeReduceIndexB;
          indexMap[ removeIndexB ] =n;
        }
      }
     
      // check if moveA has been set to n
      if( moveB < n ) {
        // check if removeIndexA has been set to n
        if( indexMap[removeIndexA] != n ) {
#ifdef DDEBUG
printf("4c\n");
#endif
          reduceIndexMap[ removeReduceIndexA ] = reduceIndexMap[ indexMap[ moveB ] ];
          indexMap[ moveB ] = removeReduceIndexA;
          indexMap[ removeIndexA ] =n;
        } else {
#ifdef DDEBUG
printf("4d\n");
#endif
          reduceIndexMap[ removeReduceIndexB ] = reduceIndexMap[ indexMap[ moveB ] ];
          indexMap[ moveB ] = removeReduceIndexB;
          indexMap[ removeIndexB ] =n;
        }
      }
    
    } 
        
    reduceIndexMap[ n - i - 1 ] =n;
    reduceIndexMap[ n - i - numberIndexToRemove ] =n;

    // accelerate termination if more than one element is terminated in a single round
    if( numberIndexToRemove == 2) i++;

  } // finish iteration
 
  /* finish */ 

#ifdef CLI 
//printf("----------------------- %d  (%d, %d)---------------\n", (int) i, (int) j, (int) sampled); 
printf(" %d:\tpoint tree\ttree\tindex\treduce\tpi\n", (int) i); 
for( l = 0; l < n; l ++) 
printf(" %d:\t%p\t%d\t%d\t%d\t%f\n", 
(int) l, 
(void *) myTree->pointerIndex[l],
(int) *( myTree->pointerIndex[l] ),
(int) ( indexMap[l] ),
(int) ( reduceIndexMap[l] ),
pi[l]
);


#ifdef DDEBUG 
printf(" --------- Print Tree -----------\n" );
printTree( myTree, myTree->root ); 
#endif

#endif

  /* delete tree */
  deleteTree(myTree);

  // free indexes */
  free(reduceIndexMap);
  free(indexMap);
  
  #ifdef CLI
    free(x);
    free(pi);
    return 0;
  #else
  
    PutRNGstate();
    return;
  #endif
}





