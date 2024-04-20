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
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

int main(void) {
  printf("\n\n    V : Mesh and displacement norm \n");
  printf("    D : Domains \n");
  printf("    N : Next domain highlighted\n\n\n");

  //
  //  -1- Lecture des donnees
  //

  femGeo *theGeometry = geoGetGeometry();
  geoMeshRead("../../data/meshSimple.txt");

  femProblem *theProblem = femElasticityRead(theGeometry, "../../data/problemSimple.txt");
  double *theSoluce = theProblem->soluce;
  // printf("System size : %d\n", theProblem->system->size);
  int n = theGeometry->theNodes->nNodes;
  femSolutiondRead(2*n,theSoluce, "../../data/UVSimple.txt");
  femElasticityPrint(theProblem);

  double *uv = malloc(2 * n * sizeof(double));
  calculateAnalytic(theProblem, uv, n, 3.e6);  // calcule la solution analytique
  for (int i = 0; i < 50; i++) {
    printf("Node %d\n", i);
    printf("Analytic solution : U = %14.7le, V = %14.7le\n", uv[2 * i], uv[2 * i + 1]);
    printf("Computed solution : U = %14.7le, V = %14.7le\n", theSoluce[2 * i], theSoluce[2 * i + 1]);
  }


  // tensions aux noeuds
  
  // X c'est R et Y c'est Z
  // calcul du tenseur des déformations infinitésimales

  double *E_XX = malloc(n * sizeof(double));
  double *E_YY = malloc(n * sizeof(double));
  double *E_XY = malloc(n * sizeof(double));
  double *E_ThetaTheta = malloc(n * sizeof(double));
  
  femElasticityEpsilon(theProblem, E_XX, 0);
  femElasticityEpsilon(theProblem, E_YY, 1);
  femElasticityEpsilon(theProblem, E_XY, 2);
  femElasticityEpsilon(theProblem, E_ThetaTheta, 3);

  // Print les epsilon
  double E = theProblem->E;
  double nu = theProblem->nu;
  double p = 3.e6;
  double E_XX_analytic = nu / E * p;
  double E_YY_analytic = -p/E;
  double E_XY_analytic = 0.0;
  double E_ThetaTheta_analytic = E_XX_analytic;
  printf(" ==== Analytic solution : E_XX = %14.7le, E_YY = %14.7le, E_XY = %14.7le, E_ThetaTheta = %14.7le\n", E_XX_analytic, E_YY_analytic, E_XY_analytic, E_ThetaTheta_analytic);
  for (int i = 0; i < 10; i++) {
    printf("Node %d\n", i);
    printf("E_XX = %14.7le, E_YY = %14.7le, E_XY = %14.7le, E_ThetaTheta = %14.7le\n", E_XX[i], E_YY[i], E_XY[i], E_ThetaTheta[i]);
    
  }
  // calcul du tenseur des contraintes

  double *sigma_XX = malloc(n * sizeof(double));
  double *sigma_ThetaTheta = malloc(n * sizeof(double));
  double *sigma_YY = malloc(n * sizeof(double));
  double *sigma_XY = malloc(n * sizeof(double));
  double *eigen_values = malloc(3 * n * sizeof(double));  // les valeurs propres
  double A = theProblem->A;
  double B = theProblem->B;
  double C = theProblem->C;
  for (int i = 0; i < n; i++) {
    sigma_XX[i] = A * E_XX[i] + B * (E_ThetaTheta[i] + E_YY[i]);
    sigma_ThetaTheta[i] =  A * E_ThetaTheta[i] + B * (E_XX[i] + E_YY[i]);
    sigma_YY[i] = A * E_YY[i] + B * (E_XX[i] + E_ThetaTheta[i]);
    sigma_XY[i] = 2 * C * E_XY[i];
  }
  // calcul des valeurs propres
  for (int i = 0; i < n; i++) {
    double *sigma = malloc(9 * sizeof(double));
    sigma[0] = sigma_XX[i];
    sigma[2] = sigma_XY[i];
    sigma[4] = sigma_ThetaTheta[i];
    sigma[6] = sigma_XY[i];
    sigma[8] = sigma_YY[i];
    double *UT = malloc(9 * sizeof(double));
    double epsilon = 1;
    int iterations = Jacobi(sigma, UT, 3, epsilon);  // l'algorithme de Jacobi, calcule les valeurs propres et les vecteurs propres
    // ordonne les valeurs propres
    if (sigma[0] > sigma[4] && sigma[0] > sigma[8]) {
      eigen_values[3 * i] = sigma[0];
      if (sigma[4] > sigma[8]) {
        eigen_values[3 * i + 1] = sigma[4];
        eigen_values[3 * i + 2] = sigma[8];
      } else {
        eigen_values[3 * i + 1] = sigma[8];
        eigen_values[3 * i + 2] = sigma[4];
      }
    } else if (sigma[4] > sigma[0] && sigma[4] > sigma[8]) {
      eigen_values[3 * i] = sigma[4];
      if (sigma[0] > sigma[8]) {
        eigen_values[3 * i + 1] = sigma[0];
        eigen_values[3 * i + 2] = sigma[8];
      } else {
        eigen_values[3 * i + 1] = sigma[8];
        eigen_values[3 * i + 2] = sigma[0];
      }
    } else {
      eigen_values[3 * i] = sigma[8];
      if (sigma[0] > sigma[4]) {
        eigen_values[3 * i + 1] = sigma[0];
        eigen_values[3 * i + 2] = sigma[4];
      } else {
        eigen_values[3 * i + 1] = sigma[4];
        eigen_values[3 * i + 2] = sigma[0];
      }
    }
    // PrintMatrix(sigma, 3);
    free(sigma);
    free(UT);
  }

  // calcul de Von Mises
  // formule générale (voir Wikipédia)
  double *sigma_VM = malloc(n * sizeof(double));  // Von Mises stress
  for (int i = 0; i < n; i++) {
    sigma_VM[i] = sqrt(0.5 * (pow(sigma_XX[i] - sigma_ThetaTheta[i], 2) + pow(sigma_ThetaTheta[i] - sigma_YY[i], 2) + pow(sigma_YY[i] - sigma_XX[i], 2) + 3 * (pow(sigma_XY[i], 2))));
  }
  double vonMisesMin = femMin(sigma_VM, n);
  double vonMisesMax = femMax(sigma_VM, n);
  printf(" ==== Minimum Von Mises stress      : %14.7e \n", vonMisesMin);
  printf(" ==== Maximum Von Mises stress      : %14.7e \n", vonMisesMax);

  // calcul de Coulomb-Mohr
  double *CM = malloc(n * sizeof(double));  // 0 ou 1 selon si la contrainte est en dessous ou au dessus de la limite de Coulomb-Mohr.
  for (int i = 0; i < n; i++) {
    // droite de coulomb mohr en ax + by + c = 0
    double a = 1.314;
    double b = 1.0;
    double c = -674073802.0;
    double xA = 0.5 * (eigen_values[3 * i] + eigen_values[3 * i + 2]);  // centre du cercle de Mohr (entre sigma1 et sigma3)
    double rMohr = 0.5 * abs((eigen_values[3 * i] - eigen_values[3 * i + 2]));  // rayon du cercle de Mohr
    double dist = abs(a * xA + c) / sqrt(a * a + b * b);  // distance entre le centre du cercle et la droite de Mohr
    if (dist < rMohr) {
      CM[i] = 1.0 * (rMohr - dist);  // satisfait pas la condition de Mohr
    } else {
      CM[i] = 0.0;  // satisfait la condition de Mohr
    }

  }
  // calcul de la moyenne des contraintes sur le top domain
  // récupère les noeuds du dit domain  => fonctionne
  femDomain *topDomain = theGeometry->theDomains[2];  // top domain
  double average_value_YY_top = 0.0;
  for (int i = 0; i < topDomain->nElem; i++) {
    int edge = topDomain->elem[i];
    int node1 = theProblem->geometry->theEdges->elem[2 * edge];
    int node2 = theProblem->geometry->theEdges->elem[2 * edge + 1];
    average_value_YY_top += sigma_YY[node1] + sigma_YY[node2];
    // printf("Edge %d : node1 = %d, node2 = %d\n", edge, node1, node2);
    // printf("sigma_ZZ_node1 = %f\n", sigma_YY[node1]);
    // printf("sigma_ZZ_node2 = %f\n", sigma_YY[node2]);
  }
  average_value_YY_top /= 2 * topDomain->nElem;
  printf(" ==== Average value of sigma_YY on the top domain       :      %f\n", average_value_YY_top);

  // calcul de la moyenne des contraintes sur bottom domain
  femDomain *bottomDomain = theGeometry->theDomains[0];  // bottom domain
  double average_value_YY_bottom = 0;
  for (int i = 0; i < bottomDomain->nElem; i++) {
    int edge = bottomDomain->elem[i];
    int node1 = theProblem->geometry->theEdges->elem[2 * edge];
    int node2 = theProblem->geometry->theEdges->elem[2 * edge + 1];
    average_value_YY_bottom += sigma_YY[node1] + sigma_YY[node2];
  }
  average_value_YY_bottom /= 2 * bottomDomain->nElem;
  printf(" ==== Average value of sigma_YY on the bottom domain    :      %f\n", average_value_YY_bottom);
  /*
  // print le déplacement sur top domain
  printf(" ==== Displacement on top domain\n");
  for (int i = 0; i < topDomain->nElem; i++) {
    int edge = topDomain->elem[i];
    int node1 = theProblem->geometry->theEdges->elem[2 * edge];
    int node2 = theProblem->geometry->theEdges->elem[2 * edge + 1];
    printf("Edge %d : node1 = %d, node2 = %d\n", edge, node1, node2);
    printf("Position node1 = %f, %f\n", theGeometry->theNodes->X[node1], theGeometry->theNodes->Y[node1]);
    printf("U_node1 = %14.7le, V_node1 = %14.7le\n", theSoluce[2 * node1], theSoluce[2 * node1 + 1]);
    printf("U_node2 = %14.7le, V_node2 = %14.7le\n", theSoluce[2 * node2], theSoluce[2 * node2 + 1]);
  }
  // print le déplacement sur bottom domain
  printf(" ==== Displacement on bottom domain\n");
  for (int i = 0; i < bottomDomain->nElem; i++) {
    int edge = bottomDomain->elem[i];
    int node1 = theProblem->geometry->theEdges->elem[2 * edge];
    int node2 = theProblem->geometry->theEdges->elem[2 * edge + 1];
    printf("Edge %d : node1 = %d, node2 = %d\n", edge, node1, node2);
    printf("Position node1 = %f, %f\n", theGeometry->theNodes->X[node1], theGeometry->theNodes->Y[node1]);
    printf("U_node1 = %14.7le, V_node1 = %14.7le\n", theSoluce[2 * node1], theSoluce[2 * node1 + 1]);
    printf("U_node2 = %14.7le, V_node2 = %14.7le\n", theSoluce[2 * node2], theSoluce[2 * node2 + 1]);
  }
  // print le déplacement sur le axis domain
  femDomain *axisDomain = theGeometry->theDomains[3];  // axis domain
  printf(" ==== Displacement on axis domain\n");
  for (int i = 0; i < axisDomain->nElem; i++) {
    int edge = axisDomain->elem[i];
    int node1 = theProblem->geometry->theEdges->elem[2 * edge];
    int node2 = theProblem->geometry->theEdges->elem[2 * edge + 1];
    printf("Edge %d : node1 = %d, node2 = %d\n", edge, node1, node2);
    printf("Position node1 = %f, %f\n", theGeometry->theNodes->X[node1], theGeometry->theNodes->Y[node1]);
    printf("U_node1 = %14.7le, V_node1 = %14.7le\n", theSoluce[2 * node1], theSoluce[2 * node1 + 1]);
    printf("U_node2 = %14.7le, V_node2 = %14.7le\n", theSoluce[2 * node2], theSoluce[2 * node2 + 1]);
  }
  */
  //
  //
  //  -2- Deformation du maillage pour le plot final
  //      Creation du champ de la norme du deplacement
  //
  femNodes *theNodes = theGeometry->theNodes;
  double deformationFactor = 1.e5;  // change le facteur de déformation pour ne pas avoir quelque chose d'absurde
  double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
  printf("X[0], Y[0] = %f, %f\n", theNodes->X[0], theNodes->Y[0]);
  printf("Deplacement en X[0], Y[0] = %f, %f\n", theSoluce[0], theSoluce[1]);
  for (int i = 0; i < n; i++) {
    normDisplacement[i] = sqrt(theSoluce[2 * i + 0] * theSoluce[2 * i + 0] + theSoluce[2 * i + 1] * theSoluce[2 * i + 1]);
  }
  for (int i = 0; i < n; i++) {
    theNodes->X[i] += deformationFactor * theSoluce[2 * i + 0];
    theNodes->Y[i] += deformationFactor * theSoluce[2 * i + 1];
  }

  double hMin = femMin(normDisplacement, n);
  double hMax = femMax(normDisplacement, n);
  printf(" ==== Minimum displacement          : %14.7e \n", hMin);
  printf(" ==== Maximum displacement          : %14.7e \n", hMax);

  //
  //  -3- Visualisation
  //
  for (int i = 0; i < n; i++) {
    sigma_YY[i] = abs(sigma_YY[i]); // mets en valeur absolue
  }

  int mode = 1;
  int domain = 0;
  int freezingButton = FALSE;
  double t, told = 0;
  char theMessage[MAXNAME];

  GLFWwindow *window = glfemInit("EPL1110 : Project 2022-23 ");
  glfwSetWindowPos(window, 0, 0);
  glfwMakeContextCurrent(window);
  glfwSetScrollCallback(window, scroll_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);

  bool autoSwtich = false;  // ce qu'il faut mettre pour avoir la "vidéo"

  if (autoSwtich) {  // ça marche

    for (mode; mode < 5; mode++) {
      double start = glfwGetTime();
      while (glfwGetTime() - start < 2.0) {
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glfemReshapeWindows(window, theGeometry->theNodes, w, h);
        t = glfwGetTime();
        if (mode == 1) {
          glfemPlotField(theGeometry->theElements, normDisplacement);
          glfemPlotMesh(theGeometry->theElements);
          sprintf(theMessage, "Deplacement elastique (* %d)  Number of elements : %d ", (int) deformationFactor, theGeometry->theElements->nElem);
          glColor3f(1.0, 0.0, 0.0);
          glfemMessage(theMessage);
        }
        if (mode == 2) // affichage de Von mises
        {
          glfemPlotField(theGeometry->theElements, sigma_VM);
          glfemPlotMesh(theGeometry->theElements);
          sprintf(theMessage, "Contrainte de Von Mises", vonMisesMax / 1.e6);
          glColor3f(1.0, 0.0, 0.0);
          glfemMessage(theMessage);
        }
        if (mode == 3) // affichage de sigma_YY (ZZ stress)
        {
          glfemPlotField(theGeometry->theElements, sigma_YY);
          glfemPlotMesh(theGeometry->theElements);
          sprintf(theMessage, "sigma_zz (moyen sur top domain : %f MPa)", average_value_YY_top / 1.e6);
          glColor3f(1.0, 0.0, 0.0);
          glfemMessage(theMessage);
          sprintf(theMessage, "sigma_zz (moyen sur bottom domain : %f MPa)", average_value_YY_bottom / 1.e6);
          glfemMessage2(theMessage);
        }
        if (mode == 4) // affichage de Coulomb-Mohr
        {
          glfemPlotField(theGeometry->theElements, CM);
          glfemPlotMesh(theGeometry->theElements);
          sprintf(theMessage, "Coulomb-Mohr respecte ou non ? Bleu : oui ; Pas bleu : non");
          glColor3f(1.0, 0.0, 0.0);
          glfemMessage(theMessage);
        }
        glfwSwapBuffers(window);
        glfwPollEvents();
      }
    }
  } else {
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
      if (glfwGetKey(window, 'Q') == GLFW_PRESS) {  // Von mises
        mode = 2;
      }
      if (glfwGetKey(window, 'W') == GLFW_PRESS) {  // ZZ stress
        for (int i = 0; i < n; i++) {
          sigma_YY[i] = -sigma_YY[i]; // mets en valeur absolue
        }
        mode = 3;
      }
      if (glfwGetKey(window, 'R') == GLFW_PRESS) {  // Coulomb-Mohr respecté ou non
        mode = 4;
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
        sprintf(theMessage, "Deplacement elastique   Number of elements : %d ", theGeometry->theElements->nElem);
        glColor3f(1.0, 0.0, 0.0);
        glfemMessage(theMessage);
      }
      if (mode == 2) // affichage de Von mises
      {
        glfemPlotField(theGeometry->theElements, sigma_VM);
        glfemPlotMesh(theGeometry->theElements);
        sprintf(theMessage, "Von mises stress", vonMisesMax / 1.e6);
        glColor3f(1.0, 0.0, 0.0);
        glfemMessage(theMessage);
      }
      if (mode == 3) // affichage de sigma_YY (ZZ stress)
      {
        glfemPlotField(theGeometry->theElements, sigma_YY);
        glfemPlotMesh(theGeometry->theElements);
        sprintf(theMessage, "sigma_zz (moyen sur top domain : %f MPa)", average_value_YY_top / 1.e6);
        glColor3f(1.0, 0.0, 0.0);
        glfemMessage(theMessage);
        sprintf(theMessage, "sigma_zz (moyen sur bottom domain : %f MPa)", average_value_YY_bottom / 1.e6);
        glfemMessage2(theMessage);
      }
      if (mode == 4) // affichage de Coulomb-Mohr
      {
        glfemPlotField(theGeometry->theElements, CM);
        glfemPlotMesh(theGeometry->theElements);
        sprintf(theMessage, "Coulomb-Mohr respecte ou non ? Bleu : oui ; Pas bleu : non");
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
  }
  

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
  free(eigen_values);
  free(sigma_VM);
  free(CM);
  femElasticityFree(theProblem);
  geoFree();
  glfwTerminate();

  exit(EXIT_SUCCESS);
  return 0;
}
