/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Preprocesseur
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "glfem.h"

int main(void) {

  //
  //  -1- Construction de la geometrie
  //
  geoInitialize();
  femGeo *theGeometry = geoGetGeometry();

  
  // OPTION 1 : Utilisation de GMSH avec OpenCascade, c'est ce que nous avons choisi d'utiliser
  theGeometry->h1 = 4.0;  // petit changement par rapport à la version donnée par les professeurs. On a 2 h donc h1 et h2 au lieu d'uniquement h
  theGeometry->h2 = 4.0;  // suite
  theGeometry->l1 = 1.0;
  theGeometry->l2 = 2.0;
  theGeometry->r = 0.5;  // < h2 et < l1


  geoMeshGenerate();
  geoMeshImport();
  
  // définition des domaines
  geoSetDomainName(16, "Bottom");
  geoSetDomainName(17, "Symmetry");
  geoSetDomainName(18, "Top");
  //Ligne supérieure précédant toute courbure 
  geoSetDomainName(8, "Upper_line_black");
  geoSetDomainName(9, "Upper_line_brown");
  //Deux courbures convexes de notre maillage 
  //Supérieure
  geoSetDomainName(2,"Upper_curvature"); 
  //Inférieure
  geoSetDomainName(5,"Lower_curvature");
  //Courbure concave de notre maillage
  geoSetDomainName(14,"Concave_curvature");
  //Partie inférieure à la courbure concave
  geoSetDomainName(10,"Lower_line_purple");
  geoSetDomainName(11,"Lower_line_brown");

  // Ecriture du maillage dans le fichier texte utilisé par PROJECT  
  geoMeshWrite("../data/mesh.txt");
  
  
  //
  //  -2- Definition du probleme
  //
  // Matériau : Carbure de Tungstène dur. Source : https://fr.wikipedia.org/wiki/Carbure_de_tungstène
  double E = 700.e9;
  double nu = 0.31;
  double rho = 15.63e3;
  double gx = 0;
  double gy = -9.81;

  femProblem *theProblem = femElasticityCreate(theGeometry, E, nu, rho, gx, gy, PLANAR_STRAIN);  // Déformation plane, pas contrainte plane
  // faut rajouter les conditions aux limites
  femElasticityAddBoundaryCondition(theProblem, "Upper_line_brown", DIRICHLET_X, 0.0, NAN);
  femElasticityAddBoundaryCondition(theProblem, "Lower_line_brown", DIRICHLET_X, 0.0, NAN);

  // femElasticityAddBoundaryCondition(theProblem, "SomethingElse", DIRICHLET_Y, 0.0, NAN);
  femElasticityPrint(theProblem);
  femElasticityWrite(theProblem, "../data/problem.txt");

  //
  //  -3- Champ de la taille de référence du maillage (uniquement pour la visualisation)
  //

  double *meshSizeField = malloc(theGeometry->theNodes->nNodes * sizeof(double));
  femNodes *theNodes = theGeometry->theNodes;
  for (int i = 0; i < theNodes->nNodes; ++i)
    meshSizeField[i] = theGeometry->geoSize(theNodes->X[i], theNodes->Y[i]);
  double hMin = femMin(meshSizeField, theNodes->nNodes);
  double hMax = femMax(meshSizeField, theNodes->nNodes);
  printf(" ==== Global requested h1 : %14.7e \n", theGeometry->h1);  // petit changement par rapport à la version donnée par les professeurs. On a 2 h donc h1 et h2 au lieu d'uniquement h
  printf(" ==== Global requested h2 : %14.7e \n", theGeometry->h2);  // suite
  printf(" ==== Minimum h          : %14.7e \n", hMin);
  printf(" ==== Maximum h          : %14.7e \n", hMax);

  //
  //  -4- Visualisation
  //

  int mode = 1;
  int domain = 0;
  int freezingButton = FALSE;
  double t, told = 0;
  char theMessage[MAXNAME];

  GLFWwindow *window = glfemInit("EPL1110 : Project 2022-23 ");
  glfwMakeContextCurrent(window);
  glfwSetScrollCallback(window, scroll_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);

  do {
    int w, h;
    glfwGetFramebufferSize(window, &w, &h);
    glfemReshapeWindows(window, theGeometry->theNodes, w, h);

    t = glfwGetTime();
    if (glfwGetKey(window, 'D') == GLFW_PRESS) {
      mode = 0;
    }
    if (glfwGetKey(window, 'V') == GLFW_PRESS) {
      mode = 1;
    }
    if (glfwGetKey(window, 'N') == GLFW_PRESS && freezingButton == FALSE) {
      domain++;
      freezingButton = TRUE;
      told = t;
    }

    if (t - told > 0.5) {
      freezingButton = FALSE;
    }
    if (mode == 1) {
      glfemPlotField(theGeometry->theElements, meshSizeField);
      glfemPlotMesh(theGeometry->theElements);
      sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
      glColor3f(1.0, 0.0, 0.0);
      glfemMessage(theMessage);
    }
    if (mode == 0) {
      domain = domain % theGeometry->nDomains;
      glfemPlotDomain(theGeometry->theDomains[domain]);
      sprintf(theMessage, "%s : %d ", theGeometry->theDomains[domain]->name, domain);
      glColor3f(1.0, 0.0, 0.0);
      glfemMessage(theMessage);
    }

    glfwSwapBuffers(window);
    glfwPollEvents();
  } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) != 1);

  // Check if the ESC key was pressed or the window was closed

  free(meshSizeField);
  femElasticityFree(theProblem);
  geoFree ();
  glfwTerminate();

  exit(EXIT_SUCCESS);
  return 0;
}
