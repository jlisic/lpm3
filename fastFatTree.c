#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>

//#include "size_tQueue.h"

/* create a node */
struct node {
  size_t dim;
  size_t * index;
  size_t indexUsed;
  double split;       /* point split on for dim */
  double * min;
  double * max;
  struct node * head; /* used for delete heavy trees */
  struct node * left;
  struct node * right; 
};

/* create typedef */
typedef struct node node;
typedef struct node * nodePtr;


/* tree node */
struct rootNode {
  size_t K;
  size_t leafSize;
  size_t n;       // rows in data
  size_t ** pointerIndex;
  double * data;
  nodePtr root;
};

/* create typedef */
typedef struct rootNode rootNode;
typedef struct rootNode * rootNodePtr;



/* function to create a new Tree */
rootNodePtr createTree( size_t K, size_t leafSize, size_t n, double * data );

/* add element to tree */
void buildTree( rootNodePtr r, nodePtr c, nodePtr p ); 

/* qsort compraison function for doubles */
int comp ( const void * a, const void * b ); 

/* create Node */
nodePtr createNode( rootNodePtr r, nodePtr p ); 

/* create children */
nodePtr * createChildren( rootNodePtr r, nodePtr c, nodePtr p); 

/* function to get the closest neighbor */
size_t getClosest( rootNodePtr r, nodePtr c, size_t item, double * dist  ); 

/* funciton to find neighbors */
void find_nn_notMe( rootNodePtr r, nodePtr c, size_t item, double * dist, size_t * query, double * queryPoint  ); 


/* return new prob */
void updateProb( double * xPtr, double * yPtr, size_t * xIndex, size_t * yIndex, size_t n ); 


/* function to create a new Tree */
rootNodePtr createTree( size_t K, size_t leafSize, size_t n, double * data ) {
  
  rootNodePtr y = malloc( sizeof( rootNode ) );
  
  y->pointerIndex = calloc( n, sizeof( size_t * ) );

  /* setup root node */
  y->K = K;
  y->leafSize = leafSize;
  y->root = NULL;
  y->data = data;
  y->n = n;

  return(y);
}




/* add an element to the tree */
void buildIndex( rootNodePtr r, nodePtr c, nodePtr p ) {
 
  size_t i; 

  /* for the root node */
  if( c == NULL) {
    c = createNode( r, NULL);
    r->root = c;
    c->indexUsed = r->n;
    for(i = 0; i < r->n; i++) c->index[i] = i; 
  }
  nodePtr * children = NULL; // place holder


  /* do we have too many points? */
  if( c->indexUsed < r->leafSize ) {

    /* save the final pointer locations */
    for( i = 0; i < c->indexUsed; i++) 
      r->pointerIndex[ c->index[i] ] = &( c->index[i] );

    return;
  } 

  /* if we have too many points 
   * create children
   * figure out our new dim
   * figure out what the split is
   * move current contents to new children
   */

//printf(" so you have decided to have children\n");
  children = createChildren( r, c, p);
//printf(" success!, they have left the house off on their own!\n");
  c->left  = children[0];
  c->right = children[1];
  free(children);

  buildIndex( r, c->left, c );  
  buildIndex( r, c->right, c );   

  return;
}




/* double comparison function for qsort */
int comp ( const void * a, const void * b ) {

  double aDbl = *((double*) a);
  double bDbl = *((double*) b);

  if( aDbl > bDbl) return 1;
  if( aDbl < bDbl) return -1;

  return 0;
}



/* function to create a simple node */
nodePtr createNode( rootNodePtr r, nodePtr p ) {

  size_t i;

  nodePtr c = malloc( sizeof( node ) );

  if( p != NULL) {
    c->index = calloc( p->indexUsed/2 +1, sizeof(size_t));
  } else {
    c->index = calloc( r->n, sizeof(size_t));
  }
  c->head  = p;
  c->left  = NULL;
  c->right = NULL;
  c->indexUsed = 0;
  c->split = 0;
  c->dim = 0;

  c->min = calloc( r->K, sizeof(double));
  c->max = calloc( r->K, sizeof(double));
  for( i = 0; i < r->K; i++) {
    (c->max)[i] = INFINITY;
    (c->min)[i] = -INFINITY;
  }
  return(c);
}




/* split and create children from a *full* node */
nodePtr * createChildren( rootNodePtr r, nodePtr c, nodePtr p) {

  double * x = NULL;
  size_t i;
  nodePtr * children = NULL;

  if( p != NULL) {
    c->dim = (p->dim + 1) % r->K; //increment mod K 
  } else {
    c->dim = 0;
  }



  /***********************************************************************************************************/

  /* get the median */
  x =  calloc( c->indexUsed, sizeof(double) ); //allocate some temporary space for finding the median

  /* use quick sort to find the median */
  for( i = 0; i < c->indexUsed; i++) {
    x[i] = r->data[ c->index[i] * r->K + c->dim];
  }  

  qsort( x, c->indexUsed, sizeof( double ), comp);  

  /* update min and max */
  c->min[c->dim] = x[0];
  c->max[c->dim] = x[c->indexUsed - 1];

  // get split
  if( c->indexUsed % 2 ) {  // is not even
    c->split = x[ c->indexUsed/2];
  } else {
    c->split = (x[ c->indexUsed / 2 ] + x[ c->indexUsed / 2 - 1 ])/2.0;
  }

  free(x); 
 
  /***********************************************************************************************************/



  children = calloc( 2, sizeof( nodePtr ) );
  
  /* allocate some memory for children */ 
  children[0] = createNode( r, c );

  /* create right node from parent */
  children[1] = createNode( r, c );


  for( i = 0; i < c->indexUsed; i++) {
    if( r->data[c->index[i] * r->K + c->dim] <= c->split ) {
     children[0]->index[children[0]->indexUsed] = c->index[i];
     children[0]->indexUsed++; 
    } else {
     children[1]->index[children[1]->indexUsed] = c->index[i];
     children[1]->indexUsed++;
    }
  }

  for( i = 0; i < r->K; i ++ ) {
    if ( i == c->dim ) {
      (children[0]->max)[i] = c->split; 
      (children[1]->min)[i] = c->split; 
    } else {
      (children[0]->max)[i] = c->max[i]; 
      (children[1]->min)[i] = c->min[i]; 
    }
    (children[1]->max)[i] = c->max[i]; 
    (children[0]->min)[i] = c->min[i]; 
  }

  free( c->index );
  c->index = NULL;

  return( children );

}




/* a function to print the tree */
void printTree( rootNodePtr r, nodePtr c ) {

  size_t i;

  printf("node: %p\n", (void *) c);
  if( c->index != NULL) {
    for( i=0; i < c->indexUsed; i++) printf("%d ", (int) c->index[i]); 
  } 
  printf("\n\tleft: %p right %p (split = %f) \n", (void *) c->left, (void*) c->right, c->split );
  printf("\n  min= ");
  for( i = 0; i < r->K; i++) printf("%f ",c->min[i] );
  printf("\n  max= ");
  for( i = 0; i < r->K; i++) printf("%f ",c->max[i] );
  printf("\n");

  if( c->left ) {
    printf("left ");
    printTree( r, c->left);
  }

  if( c->right ) {
    printf("right ");
    printTree( r, c->right);
  }

}



/* function to find the minimal Euclidian distance */
size_t getClosest( rootNodePtr r, nodePtr c, size_t item, double * dist  ) {

  size_t i,j,d;
  size_t closestIndex = r->n;

  size_t K = r->K;
  double * x = r->data;
  double * y = &(r->data[item*K]);
  double currentDist;

//printf("  getClosest: c->indexUsed = %d\n", (int) c->indexUsed );

  for( i = 0; i < c->indexUsed; i++) {

    currentDist = 0;
  
    j = c->index[i]; 

//printf("  getClosest: Checking %d against %d\n", (int) j, (int) item);

    if( j  >= r->n) continue;  // check if it's a valid index 
    if( j == item ) continue;  // don't match what we are not looking for

    for( d = 0; d < K; d++) currentDist += (x[j * K + d] - y[d]) * (x[j*K + d] - y[d]);  //calculate distance
    

    if( currentDist < *dist ) {
      *dist = currentDist; 
      closestIndex = i;
    }

  }

  if( closestIndex < r->n ) {
    return( c->index[closestIndex] );
  }

  return( r->n );
}




/* find the nearest neighbor that is not a specific index */
/* bound should be first set to the value of the node you are trying to find a neighbor for */
void find_nn_notMe( rootNodePtr r, nodePtr c, size_t item, double * dist, size_t * query, double * queryPoint  ) {

//printf("Entering %p\n", c);

  double boundDist = 0;
  size_t i;
  double distMin, distMax;
  size_t queryTmp = r->n;


  /* nothing to search for */
  if( item >= r->n ) {
//printf("nothing to do \n");
    return;
  } 
 
  /* is there anything here ? */ 
  if( c->index != NULL ) {
    queryTmp = getClosest( r, c, item, dist); 
    if( queryTmp < r->n) *query = queryTmp;

    return;
  } 

  /* move on */
  /* get bound distance */

  for(i = 0; i < r->K; i++){
      
      distMin = (queryPoint[i] - c->min[i])* (queryPoint[i] - c->min[i]);
      distMax = (queryPoint[i] - c->max[i])* (queryPoint[i] - c->max[i]);

      if( distMin < distMax ) {
        boundDist += distMin;
      } else {
        boundDist += distMax;
      }
  }  
    
//printf(" boundDist = %f ( %f )\n", boundDist , *dist );
    
 
  /* boundary distance */
  if( r->data[item * r->K + c->dim] <= c->split ) {
    if( boundDist <= *dist ) find_nn_notMe( r, c->left, item, dist, query, queryPoint );  
    if( boundDist <= *dist ) find_nn_notMe( r, c->right, item, dist, query, queryPoint );   
  } else {
    if( boundDist <= *dist ) find_nn_notMe( r, c->right, item, dist, query, queryPoint );   
    if( boundDist <= *dist ) find_nn_notMe( r, c->left, item, dist, query, queryPoint );  
  }


}


/* return new prob */
void updateProb( double * xPtr, double * yPtr, size_t * xIndex, size_t * yIndex, size_t n ) { 

  double xySum = *xPtr + *yPtr;
  double U = (double) rand()/ (double) RAND_MAX ; 

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







/* simple test case */

int main ( void )  {

  size_t n = 100000;
  size_t m = 40; // leaf size
  size_t K = 2; // dimension

  size_t i,j,l; // incrementor
  size_t k = n;
  size_t sampled;

size_t numberIndexToRemove;    // number of indexes to remove
size_t removeIndexA;           // first index to remove
size_t removeIndexB;           // second index to remove
size_t removeReduceIndexA;     // first reduced index to change 
size_t removeReduceIndexB;     // second reduced index to change 
size_t numberChosenIndexAtEnd; // number of indexes to move
size_t moveA;                  // first index to move 
size_t moveB;                  // second index to move

  double dist = INFINITY;

  /* allocate some memory and make a junk data set */ 
  double * x = calloc( n * K, sizeof(double));
  double * pi = calloc( n , sizeof(double));


  size_t * reduceIndexMap = (size_t *) calloc( n , sizeof( size_t ) );
  size_t * indexMap = (size_t *) calloc( n , sizeof( size_t ) );
  for( i=0; i< n; i++) reduceIndexMap[i]=i, indexMap[i]=i;

  //srand(0);
  srand( time(NULL) );

  /***************************** GENERATE FAKE DATA ****************************/
  for( i = 0; i < n; i ++) {
    x[i*K]    = (double) i;
    x[i*K +1] = (double) i;
    pi[i] = (double) 1000/n;
  }
  /*************************END GENERATE FAKE DATA ****************************/

  rootNodePtr myTree = createTree( K, m, n, x);

  buildIndex( myTree, NULL, NULL); 

/*`
  for( i = 0; i < n; i ++) 
    printf(" %d: %p %d %f\n", 
        (int) i, 
        (void *) myTree->pointerIndex[i],
        (int) *( myTree->pointerIndex[i] ),
        pi[i]
        );
        */
  
  /* algorithm */
  for( i = 0; i < n; i++) {

    /* randomly select j */
    sampled =  (size_t) (n - i) * ( (double) rand() / (double) RAND_MAX ) ;
    j = reduceIndexMap[sampled]; // j is the original index of the sampled data

    if( j >=n  ) {
      printf("breaking on iteration %d, for k\n", (int) i);
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

#ifdef DDEBUG
printf("----------------------- %d  (%d, %d)---------------\n", (int) i, (int) j, (int) sampled); 
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
 
 
  free(reduceIndexMap);
  free(indexMap);
  free(x);
  return 0;
}





