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
    // printf("toReturn = %f\n", toReturn);
    // toReturn = hstar;
    return fabs(toReturn);
    */
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

    // Dimensions de notre Notch, il y en aura toujours 1 
    double xNotch = (l2/2.0)+l1-r;
    double yNotch = h2 - r;
    double rNotch = r;
    //coordonnée du centre de notre 1er arc de cercle en haut
    double x_first_arc = (l2/2.0)+l1-r;
    double y_first_arc = h2+r;
    double r_first_arc = r;
    //coordonnée du centre de notre 2eme arc de cercle en bas
    double x_second_arc = l2/4.0;
    double y_second_arc = l2/4.0;
    double r_second_arc = l2/4.0;

 //
//  -1- Construction de la g�om�trie avec OpenCascade
//
 
    int ierr;
    int idPlate = gmshModelOccAddRectangle(0, 0, 0, (l2 + 2 * l1)/2.0, h1 + h2, -1, 0, &ierr);   // c'est x, y et z lower left corner puis dx, dy, tag et rounded_radius
    ErrorGmsh(ierr);
    int idRect_1 = gmshModelOccAddRectangle(l2/2.0, 0, 0, l1, h2-r, -1, 0, &ierr);   // c'est x, y et z lower left corner puis dx, dy, tag et rounded_radius
    ErrorGmsh(ierr);
    int idDiskl = gmshModelOccAddDisk(xNotch, yNotch, 0, rNotch, rNotch, -1, NULL,0,NULL,0,&ierr);  // c'est xc, yc, zc, rx, ry, tag
    ErrorGmsh(ierr);
    int idRect_2 = gmshModelOccAddRectangle(xNotch, yNotch, 0, r, r, -1, 0, &ierr);   // c'est x, y et z lower left corner puis dx, dy, tag et rounded_radius
    ErrorGmsh(ierr);
    int higherRect = gmshModelOccAddRectangle(xNotch, yNotch+r, 0, r, r, -1, 0, &ierr);   // c'est x, y et z lower left corner puis dx, dy, tag et rounded_radius
    ErrorGmsh(ierr);
    int lowerRect = gmshModelOccAddRectangle(l2/4.0, 0, 0, r_second_arc, r_second_arc, -1, 0, &ierr);   // c'est x, y et z lower left corner puis dx, dy, tag et rounded_radius
    ErrorGmsh(ierr);
    int plate[] = {2, idPlate};
    int rightBigRectangle []= {2, idRect_1};
    int Diskl[] = {2, idDiskl};
    int rightSmallRectangle[] = {2, idRect_2};
    int HigherRect[] = {2, higherRect};
    int LowerRect[] = {2, lowerRect};

    //Génération des arcs des cercles, en premier, l'arc de cercle en haut à droite
    int startTag = gmshModelOccAddPoint(x_first_arc, y_first_arc-r, 0.0, 0, -1, &ierr);
    ErrorGmsh(ierr);
    int endTag = gmshModelOccAddPoint(x_first_arc+r, y_first_arc, 0.0, 0, -1, &ierr);
    ErrorGmsh(ierr);
    int middleTag = gmshModelOccAddPoint(x_first_arc,y_first_arc, 0.0, 0, -1, &ierr);
    ErrorGmsh(ierr);
    int circleArc = gmshModelOccAddCircleArc(startTag, middleTag, endTag, -1,1, &ierr);
    ErrorGmsh(ierr);
    int line1 = gmshModelOccAddLine(endTag,middleTag, -1, &ierr);
    ErrorGmsh(ierr);
    int line2 = gmshModelOccAddLine(middleTag,startTag, -1, &ierr);
    ErrorGmsh(ierr);
    int circleArc1[] = {2, circleArc};
    int Curve_1[3] = {line1, line2, circleArc};
    int closeCurve_1 = gmshModelOccAddCurveLoop(Curve_1, 3, -1, &ierr);
    ErrorGmsh(ierr);
    int surface_1[1] = {closeCurve_1};
    int surface_1_tag = gmshModelOccAddPlaneSurface(surface_1, 1, -1, &ierr);
    ErrorGmsh(ierr);
    //Génération du deuxième arc de cercle qui est celui en bas à droite
    int startTag2 = gmshModelOccAddPoint(x_second_arc, 0, 0.0, 0, -1, &ierr);
    ErrorGmsh(ierr);
    int endTag2 = gmshModelOccAddPoint(x_second_arc+r_second_arc,y_second_arc, 0.0, 0, -1, &ierr);
    ErrorGmsh(ierr);
    int middleTag2 = gmshModelOccAddPoint(x_second_arc, y_second_arc, 0.0, 0, -1, &ierr);
    ErrorGmsh(ierr);
    int circleArc2 = gmshModelOccAddCircleArc(startTag2, middleTag2, endTag2, -1,1, &ierr);
    ErrorGmsh(ierr);
    int line3 = gmshModelOccAddLine(endTag2,middleTag2, -1, &ierr);
    ErrorGmsh(ierr);
    int line4 = gmshModelOccAddLine(middleTag2,startTag2, -1, &ierr);
    ErrorGmsh(ierr);
    int circleArc2_1[] = {2, circleArc2};
    int Curve_2[3] = {line3, line4, circleArc2};
    int closeCurve_2 = gmshModelOccAddCurveLoop(Curve_2, 3, -1, &ierr);
    ErrorGmsh(ierr);
    int surface_2[1] = {closeCurve_2};
    int surface_2_tag = gmshModelOccAddPlaneSurface(surface_2, 1, -1, &ierr);
    ErrorGmsh(ierr);
    //Création des domaines de notre projet (voir rapport pour plus de détails)
    // Ligne supérieure droite partie noire
    int StartTag_3 = gmshModelOccAddPoint((l2/2.0)+l1, h2+h1, 0.0, 0, -1, &ierr);
    int EndTag_3 = gmshModelOccAddPoint((l2/2.0)+l1, h2+h1/2.0, 0.0, 0, -1, &ierr);
    int line_4 = gmshModelOccAddLine(StartTag_3, EndTag_3, -1, &ierr);
    // Ligne supérieure droite partie brune 
    int EndTag_4 = gmshModelOccAddPoint((l2/2.0)+l1, h2+r, 0.0, 0, -1, &ierr);
    int line_5 = gmshModelOccAddLine(EndTag_3, EndTag_4, -1, &ierr);
    //Ligne inférieure droite partie mauve 
    int StartTag_5 = gmshModelOccAddPoint((l2/2.0), h2-r, 0.0, 0, -1, &ierr);
    int EndTag_5 = gmshModelOccAddPoint((l2/2.0), h2/1.5-r, 0.0, 0, -1, &ierr);
    int line_6 = gmshModelOccAddLine(StartTag_5, EndTag_5, -1, &ierr);
    //Ligne inférieure droite partie brune 
    int EndTag_6 = gmshModelOccAddPoint((l2/2.0), r, 0.0, 0, -1, &ierr);
    int line_7 = gmshModelOccAddLine(EndTag_5, EndTag_6, -1, &ierr);


    gmshModelOccCut(plate, 2, Diskl, 2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);   // c'est objectDimTags, objectDimTags_n, toolDimtags, toolDimtags_n
    ErrorGmsh(ierr);
    gmshModelOccCut(plate, 2, rightBigRectangle, 2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    ErrorGmsh(ierr);
    gmshModelOccCut(plate, 2, rightSmallRectangle, 2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    ErrorGmsh(ierr);
    gmshModelOccCut(plate, 2, HigherRect, 2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    ErrorGmsh(ierr);
    gmshModelOccCut(plate, 2, LowerRect, 2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
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
