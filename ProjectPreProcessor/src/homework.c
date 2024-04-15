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
    double y_second_arc = h2/6.3;
    double r_second_arc = r;
    //Point en plus sur le côté de notre maillage 
    double x_side = 0.0;
    double y_side = h2/5.0;
    //Paramètres de notre maille les 2 arcs ont exactement les mêmes paramètres à l'exception de leur position
    double h = ((l2/2.0)+l1)*0.3;
    double h_circle_arc = h*0.2;
    double dArc_1 = h*1.0;
    double dArc_2 = h*4.0;
    double h_Notch = h*0.05;
    double dNotch = h*4.0; 



    //Coéfficients de la fonction de la taille de la maille pour l'arc
    double a0 = -2*(h-h_circle_arc)/(dArc_1*dArc_1*dArc_1);
    double b0 = 3*(h-h_circle_arc)/(dArc_1*dArc_1);
    double c0 = 0;
    double d00 = h_circle_arc;
    //Coéfficients de la deuxieme fonction (Notch)
    double a2 = -2*(h-h_Notch)/(dNotch*dNotch*dNotch);
    double b2 = 3*(h-h_Notch)/(dNotch*dNotch);
    double c2 = 0;
    double d02 = h_Notch; 
   


    //Calcul de la distance entre le point et les formes
    double dist_arc_1 = dist(x, y, x_first_arc, y_first_arc)-r_first_arc;
    double dist_notch = dist(x, y, xNotch, yNotch)-rNotch;
    double dist_arc_2 = dist(x, y, x_second_arc, y_second_arc)-r_second_arc;
    double dist_side = dist(x, y, x_side, y_side)-r;
    /*
    if (dist_arc_1 < dArc_1 && dist_notch > dNotch) {
        return a0*dist_arc_1*dist_arc_1*dist_arc_1 + b0*dist_arc_1*dist_arc_1 + c0*dist_arc_1 + d00;
    }
    if (dist_notch < dNotch && dist_arc_1 > dArc_1) {
        return a2*dist_notch*dist_notch*dist_notch + b2*dist_notch*dist_notch + c2*dist_notch + d02;
    }
    if (dist_arc_1<dArc_1 && dist_notch<dNotch) {
        double value_1 = a0*dist_arc_1*dist_arc_1*dist_arc_1 + b0*dist_arc_1*dist_arc_1 + c0*dist_arc_1 + d00;
        double value_2 = a2*dist_notch*dist_notch*dist_notch + b2*dist_notch*dist_notch + c2*dist_notch + d02;
        return fmin(value_1, value_2);
    }
    if (dist_arc_2 < dArc_2 && dist_notch > dNotch) {
        return a0*dist_arc_2*dist_arc_2*dist_arc_2 + b0*dist_arc_2*dist_arc_2 + c0*dist_arc_2 + d00;
    }
    if (dist_notch < dNotch && dist_arc_2 > dArc_2) {
        return a2*dist_notch*dist_notch*dist_notch + b2*dist_notch*dist_notch + c2*dist_notch + d02;
    }
    if(dist_arc_2<dArc_2 && dist_notch<dNotch) {
        double value_1 = a0*dist_arc_2*dist_arc_2*dist_arc_2 + b0*dist_arc_2*dist_arc_2 + c0*dist_arc_2 + d00;
        double value_2 = a2*dist_notch*dist_notch*dist_notch + b2*dist_notch*dist_notch + c2*dist_notch + d02;
        return fmin(value_1, value_2);
    }
    return h;
    /*
    if (distl < refDist) {
        double a = 2.0 * (hstar - hd) / (refDist * refDist * refDist);
        double b = 3.0 * (hd - hstar) / (refDist * refDist);
        double d = hd;
        toReturn = a * distl * distl * distl + b * distl * distl + d;
    }
    else if (distr < refDist) {
        double a = 2.0 * (hstar - hd) / (refDist * refDist * refDist);
        double b = 3.0 * (hd - hstar) / (refDist * refDist);
        double d = hd;
        toReturn = a * distr * distr * distr + b * 
        distr * distr + d;
    }
    if (y < h2) {
        toReturn = fmin(toReturn, hd * y * y + 0.03);
    }
    */

    // printf("toReturn = %f\n", toReturn);
    // toReturn = hstar;
    return h/10.0;

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
    int Symetry_line = gmshModelOccAddLine(origine, top_left_corner,-1, &ierr);
    ErrorGmsh(ierr);
    int top_right_corner = gmshModelOccAddPoint(l1+l2/2.0, h2+h1, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int Top = gmshModelOccAddLine(top_left_corner, top_right_corner,-1, &ierr);
    ErrorGmsh(ierr);
    int low_black_domain_corner = gmshModelOccAddPoint(l1+l2/2.0, h2+h1/2.0, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int top_black_domain = gmshModelOccAddLine(top_right_corner, low_black_domain_corner,-1, &ierr);
    ErrorGmsh(ierr);
    // Start curve 1 sera la fin du domaine en mauve
    int start_curve_1 = gmshModelOccAddPoint(l1+l2/2.0, h2+r, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int purple_domain = gmshModelOccAddLine(low_black_domain_corner, start_curve_1,-1, &ierr);
    ErrorGmsh(ierr);
    int center_curve_1 = gmshModelOccAddPoint(l1+l2/2.0-r, h2+r, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int end_curve_1 = gmshModelOccAddPoint(l1+l2/2.0-r, h2, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    // le début de la courbe 2 est le même que la fin de la courbe 1
    int curve_1 = gmshModelOccAddCircleArc(start_curve_1, center_curve_1, end_curve_1,-1,1, &ierr);
    ErrorGmsh(ierr);
    int center_curve_2 = gmshModelOccAddPoint(l1+l2/2.0-r, h2-r, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int end_curve_2 = gmshModelOccAddPoint(l2/2.0, h2-r, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int curve_2 = gmshModelOccAddCircleArc(end_curve_1, center_curve_2, end_curve_2,-1,1, &ierr);
    ErrorGmsh(ierr);
    int purple_domain_2_end = gmshModelOccAddPoint(l2/2.0,h2/2.0 ,0.0,0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int purple_domain_2 = gmshModelOccAddLine(end_curve_2, purple_domain_2_end,-1, &ierr);
    ErrorGmsh(ierr);
    int brown_domain_end = gmshModelOccAddPoint(l2/2.0, l2/4.0, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int brown_domain = gmshModelOccAddLine(purple_domain_2_end, brown_domain_end,-1, &ierr);
    ErrorGmsh(ierr);
    int center_curve_3 = gmshModelOccAddPoint(l2/4.0, l2/4.0, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int end_curve_3 = gmshModelOccAddPoint(l2/4.0, 0.0, 0.0, 0.0,-1, &ierr);
    ErrorGmsh(ierr);
    int curve_3 = gmshModelOccAddCircleArc(brown_domain_end, center_curve_3, end_curve_3,-1,1, &ierr);
    ErrorGmsh(ierr);
    int bottom = gmshModelOccAddLine(end_curve_3, origine,-1, &ierr);
    ErrorGmsh(ierr);
    int curve [10] = {Symetry_line, Top, top_black_domain, purple_domain, curve_1, curve_2, purple_domain_2, brown_domain, curve_3, bottom};
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
