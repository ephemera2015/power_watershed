Licence available in the file "licence.txt"

author: Camille Couprie
(c.couprie@esiee.fr)

Softwares used :

- PINK, author Michel Couprie (m.couprie@esiee.fr)
- CSPARSE, author Tim Davis (http://www.cise.ufl.edu/research/sparse/CSparse/)

Introduction
 
 This "Powerwatershed" package provides implementation of several segmentation algorithms on 2D or 3D images.   
 
The main function for the sementation is powerwatsegm and is located
in the file powerwatsegm.c

To produce the .exe file, type "make"

Description

The algorithms are described in the following article :	
Camille Couprie, Leo Grady, Laurent Najman and Hugues Talbot , "Power Watersheds: A New Image Segmentation Framework Extending Graph Cuts, Random Walker and Optimal Spanning Forest", ICCV'09, 2009

The arguments to call the powerwatsegm.exe program are the following :
 - algorithm index : 1 to 3

  1:  Maximum Spanning Forest computed by Kruskal algorithm
  2:  Powerwatersheds (p=infinite, q=2) : Maximum Spanning Forest computed by Kruskal algorithm and Random walker on plateaus
  3:  Maximum Spanning Forest computed by Prim algorithm using Red and black trees
 
 - image name (.pgm or .ppm)
In the 3D case, the pgm format contains a pgm header
P5
3 dimensions of the image
255
followed by the image data (unsigned char).

 - seeds image (.pgm)

 For 2-labels segmentation only :
     The seeds must be a grey level pgm image, such as 
     the background is black (< 100),
     the foregroud is white (>155),
     the area to segment is grey (between 100 and 155)

 For multi-seeds segmentation :
The seeds image must be an pgm image filled with 0, and the seed
pixels must be 1, 2, ... , n(<=255) for segmentation in n labels.    

- geod :  >0 : geodesic reconstruction 0 : no geodesic reconstruction.
             
The outputs are (unless output names are specified in the options) :
     
  - "mask.pgm" the segmentation mask obtained  
  - "overlay.pgm" the segmentation overlay
  - in the case of the Powerwatershed algo, a potential map "proba.pgm"  

 Examples of use

- 2D color image segmentation :

./powerwatsegm.exe -a 2 -i images/2D/241004.ppm -m images/2D/seeds_241004_MULT.pgm 

- 3D segmentation:

./powerwatsegm.exe -a 2 -i images/3D/brainchar200.pgm -s images/3D/seeds_brain200.pgm

Visualization of the 3D pgm images :

A sotfware named 3dview is provided
Example of use :

./3dview images/3D/brainchar200.pgm

Limitations : 

The algorithms provided in this first version of the software can be
executed on 2D images and 3D images below the dimensions 255*255*200
on an ordinary GPU. More memory is needed for biger images,
different implementations of the algorithms have to be considered. 

Possible extentions :

From the current implementation that do not involve unary terms, it is
not difficult to modify the graph
to add unary terms and use the code in different energy minimization
applications. The modifications needed will concern the array of
edges index called "edges" and the neighbor functions. 

