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

  // X c'est R et Y c'est Z
  // calcul du tenseur des déformations infinitésimales

  double *E_XX = malloc(n * sizeof(double));
  double *E_YY = malloc(n * sizeof(double));
  double *E_XY = malloc(n * sizeof(double));
  double *E_ThetaTheta = malloc(n * sizeof(double));
  femElasticityAssembleElementsE_XX(theProblem, E_XX);
  femElasticityAssembleElementsE_YY(theProblem, E_YY);
  femElasticityAssembleElementsE_XY(theProblem, E_XY);
  femElasticityAssembleElementsE_ThetaTheta(theProblem, E_ThetaTheta);

  // calcul du tenseur des contraintes

  double *sigma_XX = malloc(n * sizeof(double));
  double *sigma_ThetaTheta = malloc(n * sizeof(double));
  double *sigma_YY = malloc(n * sizeof(double));
  double *sigma_XY = malloc(n * sizeof(double));
  double A = theProblem->A;
  double B = theProblem->B;
  double C = theProblem->C;
  for (int i = 0; i < n; i++) {
    sigma_XX[i] = A * E_XX[i] + B * (E_ThetaTheta[i] + E_YY[i]);
    sigma_ThetaTheta[i] =  A * E_ThetaTheta[i] + B * (E_XX[i] + E_YY[i]);
    sigma_YY[i] = A * E_YY[i] + B * (E_XX[i] + E_ThetaTheta[i]);
    sigma_XY[i] = 2 * C * E_XY[i];
  }

  // calcul de Von Mises
  // formule générale (voir Wikipédia)
  double *sigma_VM = malloc(n * sizeof(double));  // Von Mises stress
  for (int i = 0; i < n; i++) {
    sigma_VM[i] = sqrt(0.5 * (pow(sigma_XX[i] - sigma_ThetaTheta[i], 2) + pow(sigma_ThetaTheta[i] - sigma_YY[i], 2) + pow(sigma_YY[i] - sigma_XX[i], 2) + 3 * (pow(sigma_XY[i], 2))));
  }
  // Affichage de quelques valeurs pour s'assurer du bon sens de la réponse
  for (int i=0; i < 5; i++) {
    printf("E_XX[%d] = %f\n", i, E_XX[i]);
    printf("E_YY[%d] = %f\n", i, E_YY[i]);
    printf("E_XY[%d] = %f\n", i, E_XY[i]);
    printf("E_ThetaTheta[%d] = %f\n", i, E_ThetaTheta[i]);
    printf("sigma_XX[%d] = %f\n", i, sigma_XX[i]);
    printf("sigma_ThetaTheta[%d] = %f\n", i, sigma_ThetaTheta[i]);
    printf("sigma_YY[%d] = %f\n", i, sigma_YY[i]);
    printf("sigma_XY[%d] = %f\n", i, sigma_XY[i]);
    printf("sigma_VM[%d] = %f\n", i, sigma_VM[i]);
  }


  double vonMisesMin = femMin(sigma_VM, n);
  double vonMisesMax = femMax(sigma_VM, n);
  printf(" ==== Minimum Von Mises stress      : %14.7e \n", vonMisesMin);
  printf(" ==== Maximum Von Mises stress      : %14.7e \n", vonMisesMax);

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
    if (glfwGetKey(window, 'S') == GLFW_PRESS) {
      mode = 2;
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
    if (mode == 2) // affichage de Von mises
    {
      glfemPlotField(theGeometry->theElements, sigma_VM);
      glfemPlotMesh(theGeometry->theElements);
      sprintf(theMessage, "Von mises max : %f MPa", vonMisesMax / 1.e6);
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
  free(E_ThetaTheta);
  free(sigma_XX);
  free(sigma_YY);
  free(sigma_XY);
  free(sigma_ThetaTheta);
  free(sigma_VM);
  femElasticityFree(theProblem);
  geoFree();
  glfwTerminate();

  exit(EXIT_SUCCESS);
  return 0;
}
