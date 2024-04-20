//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH

#include "fem.h"


double dist(double x, double y, double x0, double y0) {  // renvoie la distance entre deux points
    return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));
}

double geoSize(double x, double y){
    return 0.04;
}



#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();
    // récupérer les valeurs de theGeometry pour les utiliser dans la fonction de création de Mesh
    double h1 = theGeometry->h1;
    double l1 = theGeometry->l1;

 //
//  -1- Construction de la g�om�trie avec OpenCascade
//
 
    int ierr;
    int rectangle = gmshModelOccAddRectangle(0.0, 0.0, 0.0, l1, h1, 0.0, -1, &ierr);
    ErrorGmsh(ierr);

//
//  -2- D�finition de la fonction callback pour la taille de r�f�rence
//      Synchronisation de OpenCascade avec gmsh
//      G�n�ration du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    geoSetSizeCallback(geoSize);
                                  
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  
       
//
//  Generation de quads :-)
//
//    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  chk(ierr);
//    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  chk(ierr);
//    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  chk(ierr);
//    gmshModelMeshGenerate(2, &ierr);  
   
 
//
//  Plot of Fltk
//
//   gmshFltkInitialize(&ierr);
//   gmshFltkRun(&ierr);  chk(ierr);
//
    
}
