/* binheap.c */

/* S. Engblom 2017-02-16 */
/* J. Cullhed 2008-08-04. */

/* Binary heap needed by the NRM and NSM. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "binheap.h"


/* This is an indexed heap. If node is the index of an element in
   "data", the place in the heap is given by INDEX[node]. The place in
   the "data" held by node n in the heap is given by INDEX2[INDEX[n]].
*/

/*----------------------------------------------------------------------*/
void initialize_heap(double *data,int *INDEX,int *INDEX2,int N)
/*** ? ***/
{
  int i;
  for (i=(N-1)>>1; i>=0; i--)
    percolate_down(i,data,INDEX,INDEX2,N);
}
/*----------------------------------------------------------------------*/
void percolate_down(int n1,double *data,int *INDEX,int *INDEX2,int N)
/*** ? ***/
{
  int child;
  int node=n1;
  double key=data[node];
  int j=INDEX[node];
  int rxn;

  while ((child = (node<<1)+1)<N) {
    if(child!=N-1 && data[child+1]<data[child]) child++;

    if(data[child]<key) {
      data[node]=data[child];
      rxn=INDEX[child];
      INDEX2[rxn]=node;
      INDEX[node]=rxn;
    }
    else
      break;

    node=child;
  }

  data[node]=key;
  INDEX[node]=j;
  INDEX2[j]=node;
}
/*----------------------------------------------------------------------*/
void percolate_up(int node,double *data,int *INDEX,int *INDEX2,int N)
/*** ? ***/
{
  int parent;
  int rxn;
  int j=INDEX[node];
  double key=data[node];

  do {
    parent=(node-1)>>1;

    if(key<data[parent]) {
      rxn=INDEX[parent];
      data[node]=data[parent];
      INDEX2[rxn]=node;
      INDEX[node]=rxn;
    }
    else
      break;

    node=parent;
  } while (parent>0);

  data[node]=key;
  INDEX[node]=j;
  INDEX2[j]=node;
}
/*----------------------------------------------------------------------*/
void update(int node,double *data,int *INDEX,int *INDEX2,int N)
/*** ? ***/
{
  int parent=(node-1)>>1;

  if(node>0 && data[node]<data[parent])
    percolate_up(node,data,INDEX,INDEX2,N);
  else
    percolate_down(node,data,INDEX,INDEX2,N);
}
/*----------------------------------------------------------------------*/
/* For debug. */
int test_heap_prty(double *data,int *INDEX,int N)
{
  int i;
  for(i=0; i<(N-1)/2; i++)
    if(data[i]>data[2*i+1] || data[i]>data[2*i+2])
      return -1;
  return 0;
}
/* /\*----------------------------------------------------------------------*\/ */
/* void print_heap(double *data,int N) */
/* {	 */
/*   int i; */
/*   printf("\n"); */
/*   for(i=0; i<(N-1)/2; i++) */
/*     printf("%.10f %.10f %.10f\n",data[i],data[2*i+1],data[2*i+2]); */
/* } */
/* /\*----------------------------------------------------------------------*\/ */
/* void test_min_prty(double *data,int N) */
/* { */
/*   int i; */
/*   double min=data[0]; */
/*   for(i=1;i<N;i++) */
/*     if(data[i]<min){ */
/*       printf("Min propery violated\n");	 */
/*       break; */
/*     } */
/* } */
/* /\*----------------------------------------------------------------------*\/ */
