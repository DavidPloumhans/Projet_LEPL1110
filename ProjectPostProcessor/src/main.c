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
  printf("\n\n    V : Mesh and displacement norm \n");
  printf("    D : Domains \n");
  printf("    N : Next domain highlighted\n\n\n");

  //
  //  -1- Lecture des donnees
  //

  femGeo *theGeometry = geoGetGeometry();
  geoMeshRead("../data/mesh.txt");

  femProblem *theProblem = femElasticityRead(theGeometry, "../data/problem.txt");
  double *theSoluce = theProblem->soluce;
  int n = theGeometry->theNodes->nNodes;
  femSolutiondRead(2*n,theSoluce, "../data/UV.txt");
  femElasticityPrint(theProblem);

  double *E_XX = malloc(n * sizeof(double));
  double *E_YY = malloc(n * sizeof(double));
  double *E_XY = malloc(n * sizeof(double));
  femElasticityAssembleElementsE_XX(theProblem, E_XX);
  femElasticityAssembleElementsE_YY(theProblem, E_YY);
  femElasticityAssembleElementsE_XY(theProblem, E_XY);

  // pas sur du calcul des contraintes
  double *sigma_XX = malloc(n * sizeof(double));
  double *sigma_YY = malloc(n * sizeof(double));
  double *sigma_XY = malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) {
    sigma_XX[i] = theProblem->A * E_XX[i] + theProblem->B * E_YY[i];
    sigma_YY[i] = theProblem->A * E_YY[i] + theProblem->B * E_XX[i];
    sigma_XY[i] = 2*theProblem->C * E_XY[i];
  }
  double *sigma_VM = malloc(n * sizeof(double));  // Von Mises stress
  //
  //  -2- Deformation du maillage pour le plot final
  //      Creation du champ de la norme du deplacement
  //

  femNodes *theNodes = theGeometry->theNodes;
  double deformationFactor = 1e1;  // change le facteur de déformation pour ne pas avoir quelque chose d'absurde
  double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));

  for (int i = 0; i < n; i++) {
    theNodes->X[i] += theSoluce[2 * i + 0] * deformationFactor;
    theNodes->Y[i] += theSoluce[2 * i + 1] * deformationFactor;
    normDisplacement[i] = sqrt(theSoluce[2 * i + 0] * theSoluce[2 * i + 0] + theSoluce[2 * i + 1] * theSoluce[2 * i + 1]);
  }

  double hMin = femMin(normDisplacement, n);
  double hMax = femMax(normDisplacement, n);
  printf(" ==== Minimum displacement          : %14.7e \n", hMin);
  printf(" ==== Maximum displacement          : %14.7e \n", hMax);

  //
  //  -3- Visualisation
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
      glfemPlotField(theGeometry->theElements, normDisplacement);
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

  free(normDisplacement);
  free(E_XX);
  free(E_YY);
  free(E_XY);
  free(sigma_XX);
  free(sigma_YY);
  free(sigma_XY);
  femElasticityFree(theProblem);
  geoFree();
  glfwTerminate();

  exit(EXIT_SUCCESS);
  return 0;
}
