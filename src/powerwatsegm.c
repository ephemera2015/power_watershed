/*Copyright ESIEE (2009) 

Author :
Camille Couprie (c.couprie@esiee.fr)

Contributors : 
Hugues Talbot (h.talbot@esiee.fr)
Leo Grady (leo.grady@siemens.com)
Laurent Najman (l.najman@esiee.fr)

This software contains some image processing algorithms whose purpose is to be
used primarily for research.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/


#include <stdio.h>
#include <sys/types.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <mccodimage.h>
#include <mcimage.h>
#include <cccodimage.h>
#include <lMSF.h>
#include <MSF_RW.h>
#include <powerwatsegm.h>
#include <image_toolbox.h>
#include <unistd.h>
#include <argv.h>
#include <iostream>

/* =============================================================== */
int main(int argc, char **argv) 
/* =============================================================== */
{
  int algo =-1;
  char * image_name = NULL;
  char * seeds_name =NULL;
  char * multi_seeds_name =NULL;
  char * overlay_name =NULL;
  char * output_name =NULL;
  bool geod = false;
  argv_t args[] = {
    { 'a', 0L, ARGV_INT, &algo, (char*)" 1|2|3", (char*)"algo: Kruskal(1) PW q=2(2) Prim(3)" },
    { 'i', 0L, ARGV_CHAR_P, &image_name, (char*)" *.pgm |*.ppm",  (char*)"image_name" },
    { 's', 0L, ARGV_CHAR_P, &seeds_name, (char*)" *.pgm",  (char*)"seed_name (see README)" },
    { 'm', 0L, ARGV_CHAR_P, &multi_seeds_name, (char*)" *.pgm",  (char*)"multi_seed_name (see README)" },
    { 'g', 0L, ARGV_BOOL, &geod, (char*)" 0|1", (char*)"geod: reconstruction of the weights(1)" },
    { 'o', 0L, ARGV_CHAR_P, &output_name, (char*)" *.pgm",  (char*)"output mask name" },
    { 'v', 0L, ARGV_CHAR_P, &overlay_name, (char*)" *.ppm",  (char*)"output overlay name" },
    { ARGV_LAST, 0L, 0, 0L, 0L, 0L}
  };
 
  argv_process(args, argc, argv);
  if(seeds_name ==NULL) 
  {
    seeds_name= multi_seeds_name;
  }

  if ((algo==-1) ||(image_name ==NULL  )||(seeds_name ==NULL))
  {
      fprintf(stderr, "usage: %s -a algo -i image(.pgm or .ppm) <-s seeds_2_labels.pgm | -m seeds_m_labels.pgm > \n", argv[0]);
      fprintf(stderr, "options : [-g geod] [-o output mask name] [-v image overlay name]\n");
      fprintf(stderr, "algo : 1 : Kruskal \n");
      fprintf(stderr, "       2 : PW q=2 \n");
      fprintf(stderr, "       3 : Prim \n");
      fprintf(stderr, "seeds[MULT].pgm  (see the README for more precision) \n");
      fprintf(stderr, "geod : 1 : geodesic reconstruction \n");
      fprintf(stderr, "       0 : no reconstruction \n");
      exit(1);
  }
  
  int t1=clock();	
  bool quicksort = false;
  int32_t i=0,j=0;
  struct xvimage * image=readimage(image_name);
  struct xvimage * output = NULL;
  struct xvimage * seeds=readimage(seeds_name);

  
  int rs=rowsize(seeds), cs=colsize(seeds), ds=depth(seeds);  
  size_t  N = rs * cs * ds,M = ds*rs*(cs-1)+ds*(rs-1)*cs+(ds-1)*cs*rs;
  unsigned char * s=UCHARDATA(seeds),*img= UCHARDATA(image);
  int sum=0;
  for(i=0;i<N;++i)
  {
    if(img[i]==0)
    {
        s[i]=0;
        ++sum;
    }
  }

  int ** edges;     
  int * index_seeds = (int*)malloc(N*sizeof(int));
  uint8_t * index_labels = (uint8_t*)malloc(N*sizeof(uint8_t));
  //multilabel seed image

  int nblabels = 0;
  for (i=0;i<N;i++)
  {
    if(s[i]>0)
    {
	  index_seeds[j]=i;
	  index_labels[j]=s[i];
	  j++;
	  if(s[i]>nblabels) nblabels = s[i];
    }
  }

  int size_seeds = j;//number of seeds
  std::cout<<N-size_seeds-sum<<std::endl;
  edges =  (int**)malloc(2*sizeof(int*));
  for(i=0;i<2;i++) edges[i] = (int*)malloc(M*sizeof(int));
  compute_edges(edges,rs, cs, ds);
  if (algo == 1) // Kruskal
  {
    uint32_t * weights = (uint32_t *)malloc(sizeof(uint32_t)*M);
    int max_weight = 255;
    int cnt=grey_weights(image, weights, edges,index_seeds, size_seeds, geod, quicksort);
    output = MSF_Kruskal(img,edges,weights, max_weight, index_seeds, index_labels, size_seeds, rs, cs, ds, nblabels,sum);
    free(weights);
  }

 /*
  else if (algo == 3) // Prim RB tree
    {  
      uint32_t * weights = (uint32_t *)malloc(sizeof(uint32_t)*M);
	grey_weights(image_name, weights, edges,index_seeds, size_seeds, geod, quicksort);
      output = MSF_Prim(edges,weights, index_seeds, index_labels, size_seeds,rs, cs, ds, nblabels);
      free(weights);
    }*/

 else if (algo == 2) // Kruskal & RW on plateaus multiseeds linear time
 {    
      uint32_t * weights = (uint32_t *)malloc(sizeof(uint32_t)*M);
      uint32_t * normal_weights ;
      uint32_t max_weight = 255;
      normal_weights = grey_weights_PW(image_name,  edges,index_seeds, size_seeds, weights, quicksort);
      if (geod ==true)
      {
        //output = PowerWatershed_q2(edges, weights, weights, max_weight,index_seeds, index_labels, size_seeds,rs, cs, ds, nblabels, quicksort, img_proba);
      }
      else
      {
        output = PowerWatershed_q2(edges, weights, normal_weights,max_weight,index_seeds, index_labels, size_seeds,rs, cs, ds, nblabels, quicksort);
      }
      free(weights);
      free(normal_weights);
     
  }

  int t2=clock();
  assert(output != NULL);
 
  free(index_seeds);
  free(index_labels);
  if (output_name == NULL) 
      output_name =(char*)"mask.pgm"; 
  writeimage(output, output_name);

  // overlay for 2D images onl
  
  freeimage(image);
  freeimage(output);
  freeimage(seeds);
  printf("Computation time : %.6lf seconds elapsed\n", ((double)t2-t1)/CLOCKS_PER_SEC);
  return 0;
} 


