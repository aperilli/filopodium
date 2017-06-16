#include <stdio.h>
#include "sprng.h"  /* SPRNG header file */
#include <mpi.h>
#include <stdlib.h>

int *stream, *cm_stream;

void easy_init_rng_(int *SEED , int *myrank, int *nproc, int *idebug)
{
/* extern void *bcast_int(); */
  if ( *SEED == 0 ) 
    *SEED = make_sprng_seed();
  stream = init_sprng(*myrank,*nproc+1,*SEED,SPRNG_DEFAULT); 




          if(*idebug!=0) printf("easy_init_rng: seed before bcast= %d %d\n",*SEED,*myrank);
/*  bcast_int(*SEED,1,0); */
MPI_Bcast((void*)SEED,1,MPI_INTEGER,0,MPI_COMM_WORLD);
          if(*idebug!=0) printf("easy_init_rng: seed after bcast= %d %d\n",*SEED,*myrank);
  cm_stream = init_sprng(*nproc,*nproc+1,*SEED,SPRNG_DEFAULT);

          if (*idebug!=0) printf("first Random#: %f %d\n",sprng(stream),*myrank);
          if (*idebug!=0) printf("first CM_Random#: %f %d\n",sprng(cm_stream),*myrank);
}

void easy_print_rng_()
{
  printf(" Print information about local stream:\n");
  print_sprng(stream);
}

void easy_print_rng_cm_()
{
  printf(" Print information about common stream:\n");
  print_sprng(cm_stream);
}

double rng_()
{
  return sprng(stream); 
}

double rng_cm_()
{
  return sprng(cm_stream);
}

int irng_()
{
  return isprng(stream);
}

void pack_rng_stream_(int *size, char *pack)
{
  char *bytes;
  int size_local;
  int i;

  size_local = pack_sprng(stream,&bytes);
  *size = size_local;
  for(i=0; i<size_local; i++)
    pack[i] = bytes[i];

  free(bytes);
}

void pack_rng_cm_stream_(int *size, char *pack)
{
  char *bytes;
  int size_local;
  int i;

  size_local = pack_sprng(cm_stream,&bytes);
  *size = size_local;
  for(i=0; i<size_local; i++)
    pack[i] = bytes[i];

  free(bytes);
}

void unpack_rng_stream_(char *pack)
{
  stream = unpack_sprng(pack);
}

void unpack_rng_cm_stream_(char *pack)
{
  cm_stream = unpack_sprng(pack);
}

int get_rng_pack_length_()
{
  return MAX_PACKED_LENGTH;
}
