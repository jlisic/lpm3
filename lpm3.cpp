#include <flann/flann.hpp>
#include <flann/io/hdf5.h>
#include <stdio.h> 
#include <stdlib.h>

#include "R.h"
#include "Rinternals.h"
#define CSTACK_DEFNS 7
#include "Rinterface.h"

#include "Rmath.h"
#include "omp.h"



// function to update the probability
// Note that this is destructive with respect to prob
void updateProb( double * prob, size_t i, size_t j) {

  double probI = prob[i];
  double probJ = prob[j];
  double totalProb = probI + probJ;

  // handle the less than 1 total prob case 
  if( totalProb < 1 )  {
    if ( runif(0.0, 1.0) < probJ / totalProb )  {
      prob[i] = 0;
      prob[j] = totalProb;
      return;
    }
    prob[i] = totalProb;
    prob[j] = 0;
    return;
  }

  if ( runif(0.0, 1.0) < ( 1 - probJ) / ( 2 - totalProb ) ) {
    prob[i] = 1;
    prob[j] = totalProb -1;
    return;
  }
    
  prob[i] = totalProb -1;
  prob[j] = 1;
  return;

} 




// extern start for .C interface
extern "C" {


/* local pivotal method version 3*/
void R_lpm3(
  double * x,
  double * prob,
  int * xRowPtr,
  int * xColPtr
) {

  /* simplify the interface by dereferencing the pointers */
  size_t xRow =     (size_t) * xRowPtr; 
  size_t xCol =     (size_t) * xColPtr;
  size_t i,j,k,l;
  size_t sampled;

  /* create some space for distances */
  double * distances = (double *) calloc( 2, sizeof( double ) );
  int * neighbors = (int *) calloc( 2 , sizeof( int ) );
  
  int * reduceIndexMap = (int *) calloc( xRow , sizeof( int ) );
  int * indexMap = (int *) calloc( xRow , sizeof( int ) );
  for( l =0; l < xRow; l++) reduceIndexMap[l]=l, indexMap[l]=l;

  /* fixes for R with OpenMP */
  R_CStackLimit=(uintptr_t)-1;

  flann::Matrix<double> xMatrix;

  flann::Matrix<double> yMatrix(x, xRow, xCol); 
  flann::Matrix<double> distancesMatrix( distances, 1, 2 ); 
  flann::Matrix<int>    neighborsMatrix( neighbors, 1, 2 ); 
  
  /* get RNG State */
  GetRNGstate();


  /******************* CREATE INDEX *********************/ 
  /*
  AutotunedIndexParams:
    float target_precision = 1.0,  // exact match
    float build_weight = 0.01,
    float memory_weight = 0,
    float sample_fraction = 0.1
  */

  //flann::AutotunedIndexParams lpmIndexParameters( 1.00, 0.01, 0, 0.1);
  flann::KDTreeIndexParams lpmIndexParameters( 1 );
  flann::Index< flann::L2<double> > index( yMatrix, lpmIndexParameters  ); 

  /* build the index */ 
  index.buildIndex();

  /******************* SEARCH *********************/ 
  flann::SearchParams mpSearchParameters( -1, 0, true);
  
  //for( i = 0; i < *iterations; i++ ) {
  for( i = 0; i < xRow; i++ ) {
  //Rprintf("Iteration: %d ----------------------------------------------------------- %d \n", (int) i, (int) xRow - i - 1);
  //Rprintf("Reduce Index Map / Prob:\n");
  //for(l=0; l < xRow; l++) printf("[%d]\t %d\t %d\t  %f\n", (int) l, reduceIndexMap[l], indexMap[l] , prob[l]);
 
    // sample a value  
    sampled = (size_t) floor( runif(0.0,1.0) * (xRow - i) ) ; 
    

    // the reduce mapping maintains a listing of non-zero, non-one elements
    j = reduceIndexMap[sampled];

    xMatrix = flann::Matrix<double>(&x[j*xCol], 1, xCol); 

    /* look for points in the training set */
    index.knnSearch(xMatrix, neighborsMatrix, distancesMatrix, 2, mpSearchParameters); 

    k = neighbors[1];
    
    //Rprintf("sampled = %d, j = %d k = %d\n", (int) sampled, (int) j, (int) k  );

    // update prob
    updateProb(prob,j,k);


    // update tree and reduce mapping
    //
    // Example
    // Let 3 be selected and 5 gets the weight for both (but is still less than 0 )
    // 0 1 2 [3] 4 [5] 6
    //
    // reduceIndexMap
    // 0 1 2 6 3 5 -1
    //
    // indexMap
    // 0 1 2 -1 4 5 3 
    //
    if( prob[j] <= 0) { 
      index.removePoint(j); // remove point from tree
      if( j != xRow - i - 1 ) {
        reduceIndexMap[ sampled ] = reduceIndexMap[xRow -i - 1]; // replace j with last point
        indexMap[ reduceIndexMap[xRow -i -1] ] = sampled; 
      }
      reduceIndexMap[xRow -i - 1] = -1;
      indexMap[j] = -1;
    }
    if( prob[k] <= 0) { 
      index.removePoint(k); // remove point from tree
      if( k != xRow - i - 1 ) {
        reduceIndexMap[ indexMap[k] ] = reduceIndexMap[xRow -i - 1]; // replace k with last point
        indexMap[ reduceIndexMap[xRow -i -1] ] = indexMap[k]; 
      }
      reduceIndexMap[xRow -i - 1] = -1;
      indexMap[k] = -1;
    }

    if( prob[j] >= 1) { 
      index.removePoint(j); // remove point from tree
      if( j != xRow - i - 1 ) {
        reduceIndexMap[ sampled ] = reduceIndexMap[xRow -i - 1]; // replace j with last point
        indexMap[ reduceIndexMap[xRow -i -1] ] = sampled; 
      }
      reduceIndexMap[xRow -i - 1] = -1;
      indexMap[j] = -1;
      continue;
    }
    

    if( prob[k] >= 1) { 
      index.removePoint(k); // remove point from tree
      if( k != xRow - i - 1 ) {
        reduceIndexMap[ indexMap[k] ] = reduceIndexMap[xRow -i - 1]; // replace k with last point
        indexMap[ reduceIndexMap[xRow -i -1] ] = indexMap[k]; 
      }
      reduceIndexMap[xRow -i - 1] = -1;
      indexMap[k] = -1;
      continue;
    }

  }
 
  free( distances );
  free( neighbors );
  free(reduceIndexMap);
  free(indexMap);

  /* return the RNG State */
  PutRNGstate();

  return; 
}


// extern end for .C interface
}



