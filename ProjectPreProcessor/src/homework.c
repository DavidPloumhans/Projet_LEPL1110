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

    femGeo* theGeometry = geoGetGeometry();
    // faut faire un truc un peu opti
    //Paramètres de notre maillage 
    double h1 = theGeometry->h1;
    double h2 = theGeometry->h2;
    double l1 = theGeometry->l1;
    double l2 = theGeometry->l2;
    double r = theGeometry->r;
    double pi = 3.14159265358979323846;
    double pi_4 = pi/4.0;
    //Coordonnées de nos formes importantes
    // Dimensions de notre Notch, il y en aura toujours 1 
    double xNotch = (l2/2.0)+l1-r;
    double yNotch = h2 - r;
    double rNotch = r;
    //coordonnée du point le plus extrême de notre 1er arc de cercle en haut
    double x_first_arc = ((l2/2.0)+l1-r);//+r*cos(pi_4);
    double y_first_arc = (h2+r);//-r*sin(pi_4);
    double r_first_arc = r;
    //coordonnée du centre de notre 2eme arc de cercle en bas
    double x_second_arc = (l2/2.0)-r;
    double y_second_arc = l2/4.0;
    double r_second_arc = l2/4.0;
    //Point en plus sur le côté de notre maillage 
    double x_side = 0.0;
    double y_side = h2/5.0;
    //Paramètres de notre maille les 2 arcs ont exactement les mêmes paramètres à l'exception de leur position
    double h = ((l2/2.0)+l1)*0.3;
    double h_circle_arc = h*0.2;
    double h_Notch = h*0.05;
    double dNotch = h*3.0; 
    //Coéfficients de la deuxieme fonction (Notch)
    double a2 = -2*(h-h_Notch)/(dNotch*dNotch*dNotch);
    double b2 = 3*(h-h_Notch)/(dNotch*dNotch);
    double c2 = 0;
    double d02 = h_Notch; 
   

    //Calcul de la distance entre le point et les formes
    double dist_notch = dist(x, y, xNotch, yNotch)-rNotch;
    double to_return = h*0.45;
    if (dist_notch < dNotch) {
        to_return = fmin((a2*dist_notch*dist_notch*dist_notch + b2*dist_notch*dist_notch + c2*dist_notch + d02)*1.5, to_return);
    }

    if (y<h2) {
        to_return = fmin(to_return, (h * 0.20 * (0.3*y+0.3))*1.5);
    }
return to_return;
}



#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();
    // récupérer les valeurs de theGeometry pour les utiliser dans la fonction de création de Mesh
    double h1 = theGeometry->h1;
    double h2 = theGeometry->h2;
    double l1 = theGeometry->l1;
    double l2 = theGeometry->l2;
    double r = theGeometry->r;


 //
//  -1- Construction de la g�om�trie avec OpenCascade
//
 
    int ierr;
    int origine = gmshModelOccAddPoint(0.0, 0.0, 0.0, 0.0, -1, &ierr);  
    ErrorGmsh(ierr);
    int top_left_corner = gmshModelOccAddPoint(0.0, h2+h1, 0.0, 0.0,-1,&ierr);
    ErrorGmsh(ierr);
    int Symetry_line = gmshModelOccAddLine(top_left_corner, origine,-1, &ierr);
    ErrorGmsh(ierr);
    int top_right_corner = gmshModelOccAddPoint(l1+l2/2.0, h2+h1, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int Top = gmshModelOccAddLine(top_right_corner, top_left_corner,-1, &ierr);
    ErrorGmsh(ierr);
    int low_black_domain_corner = gmshModelOccAddPoint(l1+l2/2.0, h2+h1/2.0, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int top_black_domain = gmshModelOccAddLine(low_black_domain_corner, top_right_corner,-1, &ierr);
    ErrorGmsh(ierr);
    // Start curve 1 sera la fin du domaine en mauve
    int start_curve_1 = gmshModelOccAddPoint(l1+l2/2.0, h2+r, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int purple_domain = gmshModelOccAddLine(start_curve_1, low_black_domain_corner,-1, &ierr);
    ErrorGmsh(ierr);
    int center_curve_1 = gmshModelOccAddPoint(l1+l2/2.0-r, h2+r, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int end_curve_1 = gmshModelOccAddPoint(l1+l2/2.0-r, h2, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    // le début de la courbe 2 est le même que la fin de la courbe 1
    int curve_1 = gmshModelOccAddCircleArc(end_curve_1, center_curve_1, start_curve_1,-1,1, &ierr);
    ErrorGmsh(ierr);
    int center_curve_2 = gmshModelOccAddPoint(l1+l2/2.0-r, h2-r, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int end_curve_2 = gmshModelOccAddPoint(l2/2.0, h2-r, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int curve_2 = gmshModelOccAddCircleArc(end_curve_2, center_curve_2, end_curve_1,-1,1, &ierr);
    ErrorGmsh(ierr);
    int purple_domain_2_end = gmshModelOccAddPoint(l2/2.0,h2/2.0 ,0.0,0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int purple_domain_2 = gmshModelOccAddLine(purple_domain_2_end, end_curve_2,-1, &ierr);
    ErrorGmsh(ierr);
    int brown_domain_end = gmshModelOccAddPoint(l2/2.0, l2/4.0, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int brown_domain = gmshModelOccAddLine(brown_domain_end, purple_domain_2_end,-1, &ierr);
    ErrorGmsh(ierr);
    int center_curve_3 = gmshModelOccAddPoint(l2/4.0, l2/4.0, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int end_curve_3 = gmshModelOccAddPoint(l2/4.0, 0.0, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int curve_3 = gmshModelOccAddCircleArc(end_curve_3, center_curve_3, brown_domain_end,-1,1, &ierr);
    ErrorGmsh(ierr);
    int bottom = gmshModelOccAddLine(origine, end_curve_3,-1, &ierr);
    ErrorGmsh(ierr);
    int curve [10] = {bottom, curve_3, brown_domain, purple_domain_2, curve_2, curve_1, purple_domain, top_black_domain, Top, Symetry_line};
    int domain = gmshModelOccAddWire(curve, 10, -1,1, &ierr);
    ErrorGmsh(ierr);
    int domain_pointer[1]={domain};
    int surface = gmshModelOccAddPlaneSurface(domain_pointer,1,-1, &ierr);
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
