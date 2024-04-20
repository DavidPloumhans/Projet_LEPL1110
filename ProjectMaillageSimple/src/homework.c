#include "fem.h"
// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)

void femElasticityAssembleElements(femProblem *theProblem) {

  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->rule;
  femDiscrete *theSpace = theProblem->space;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theMesh = theGeometry->theElements;
  double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
  int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];
  int nLocal = theMesh->nLocalNode;
  double a = theProblem->A;
  double b = theProblem->B;
  double c = theProblem->C;
  double rho = theProblem->rho;
  double gx = theProblem->gx;
  double gy = theProblem->gy;
  double **A = theSystem->A;
  double *B = theSystem->B;
  int count = 0;
  for (iElem = 0; iElem < theMesh->nElem; iElem++) {
    for (j = 0; j < nLocal; j++) {
      map[j] = theMesh->elem[iElem * nLocal + j];
      mapX[j] = 2 * map[j];
      mapY[j] = 2 * map[j] + 1;
      x[j] = theNodes->X[map[j]];
      y[j] = theNodes->Y[map[j]];
    }
    
    /*
    printf("Element %d\n", iElem);
    printf("X[0] = %f, Y[0] = %f\n", x[0], y[0]);
    printf("X[1] = %f, Y[1] = %f\n", x[1], y[1]);
    printf("X[2] = %f, Y[2] = %f\n", x[2], y[2]);
    */
   
    for (iInteg = 0; iInteg < theRule->n; iInteg++) {
      
      double xsi = theRule->xsi[iInteg];
      double eta = theRule->eta[iInteg];
      double weight = theRule->weight[iInteg];
      femDiscretePhi2(theSpace, xsi, eta, phi);
      femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

      double dxdxsi = 0.0;
      double dxdeta = 0.0;
      double dydxsi = 0.0;
      double dydeta = 0.0;
      double xLoc = 0.0;
      for (i = 0; i < theSpace->n; i++) {
        dxdxsi += x[i] * dphidxsi[i];
        dxdeta += x[i] * dphideta[i];
        dydxsi += y[i] * dphidxsi[i];
        dydeta += y[i] * dphideta[i];
        xLoc += x[i]*phi[i];
      }
      double jac = dxdxsi * dydeta - dxdeta * dydxsi;
      if (jac < 0.0 && iInteg == 0) {
        printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n");
        count++;
      }
      jac = fabs(jac);

      for (i = 0; i < theSpace->n; i++) {  // Calcul des dérivées des fonctions de forme
          dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
          dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
      }
      
      if (theProblem->planarStrainStress != AXISYM) {
        // Cas non-axisymétrique
        for (i = 0; i < theSpace->n; i++) {
          for (j = 0; j < theSpace->n; j++) {
            A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight;
            A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight;
            A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight;
            A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight;
          }
        }
        for (i = 0; i < theSpace->n; i++) {
          B[mapX[i]] += phi[i] * gx * rho * jac * weight;
          B[mapY[i]] += phi[i] * gy * rho * jac * weight;
        }
      } else if (theProblem->planarStrainStress == AXISYM) {
      // Cas axisymétrique
      for (i = 0; i < theSpace->n; i++) {
        for(j = 0; j < theSpace->n; j++) {
          A[mapX[i]][mapX[j]] += (dphidx[i] * a * xLoc * dphidx[j] + dphidy[i] * c * xLoc * dphidy[j] + phi[i] * ((b * dphidx[j]) + (a * phi[j] / xLoc)) + dphidx[i] * b * phi[j]) * jac * weight;
          A[mapX[i]][mapY[j]] += (dphidx[i] * b * xLoc * dphidy[j] + dphidy[i] * c * xLoc * dphidx[j] + phi[i] * b * dphidy[j]) * jac * weight;
          A[mapY[i]][mapX[j]] += (dphidy[i] * b * xLoc * dphidx[j] + dphidx[i] * c * xLoc * dphidy[j] + dphidy[i] * b * phi[j]) * jac * weight;
          A[mapY[i]][mapY[j]] += (dphidy[i] * a * xLoc * dphidy[j] + dphidx[i] * c * xLoc * dphidx[j]) * jac * weight;
        }
      }
      for (i = 0; i < theSpace->n; i++) {         
        B[mapX[i]] += phi[i] * gx * rho * jac * weight * xLoc;
        B[mapY[i]] += phi[i] * gy * rho * jac * weight * xLoc;
      }
      } else {Error("Invalid PlanarStrainStress value. Must be PLANE_STRESS or PLANE_STRAIN or AXISYM\n");}
    }
  }
}

void femElasticityAssembleNeumann(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->ruleEdge;
  femDiscrete *theSpace = theProblem->spaceEdge;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theEdges = theGeometry->theEdges;
  double x[2], y[2], phi[2];
  int iBnd, iElem, iInteg, iEdge, i, j, d, map[2];
  int nLocal = 2;
  double *B = theSystem->B;

  for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {
    femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
    femBoundaryType type = theCondition->type;
    double value = theCondition->value1;

    if(type != NEUMANN_X && type != NEUMANN_Y && type != NEUMANN_N && type != NEUMANN_T){
      continue;
    }

    for (iEdge = 0; iEdge < theCondition->domain->nElem; iEdge++) {
      iElem = theCondition->domain->elem[iEdge];
      for (j = 0; j < nLocal; j++) {
        map[j] = theEdges->elem[iElem * nLocal + j];
        x[j] = theNodes->X[map[j]];
        y[j] = theNodes->Y[map[j]];
      }

      double tx = x[1] - x[0];
      double ty = y[1] - y[0];
      double length = hypot(tx, ty);
      double jac = length / 2.0;
      
      double f_x = 0.0;
      double f_y = 0.0;
      if (type == NEUMANN_X) {  // f_y sera nul
        f_x = value;
      }
      if (type == NEUMANN_Y) {  // f_x sera nul
        f_y = value;
      }
      // Attention, pour le normal tangent on calcule la normale (sortante) au SEGMENT, surtout PAS celle de constrainedNodes
      // Une petite aide pour le calcul de la normale :-)
      // le triangle trigonométrique est donné dans le sens antihorlogique (d'où le - pour le ny) et la normale est sortante. La tangente est dans le sens anti-horlogique.
      double nx =  ty / length;  // sinus de l'angle en x
      double ny = -tx / length;  // - cosinus de l'angle en x
      if (type == NEUMANN_N) {
        f_x = value * nx;
        f_y = value * ny;
      } else if (type == NEUMANN_T) {
        f_x = -value * ny;
        f_y = value * nx;
      }

      if (theProblem->planarStrainStress != AXISYM) {
        for (iInteg = 0; iInteg < theRule->n; iInteg++) {
          double xsi = theRule->xsi[iInteg];
          double weight = theRule->weight[iInteg];
          femDiscretePhi(theSpace, xsi, phi);
          for (i = 0; i < theSpace->n; i++) {
            B[2*map[i] + 0] += jac * weight * phi[i] * f_x;
            B[2*map[i] + 1] += jac * weight * phi[i] * f_y;
          }
        }
      } else if (theProblem->planarStrainStress == AXISYM) {  // même chose mais facteur xLoc en plus
        for (iInteg = 0; iInteg < theRule->n; iInteg++) {
          double xLoc = 0.0;  // c'était un int ici avant !!
          double xsi = theRule->xsi[iInteg];
          double weight = theRule->weight[iInteg];
          femDiscretePhi(theSpace, xsi, phi);
          // printf("theSpace->n : %d\n", theSpace->n);
          for (i = 0; i < theSpace->n; i++) {
            xLoc += x[i]*phi[i];  // x local
            // printf("x[i] : %f\n", x[i]);
            // printf("phi[i] : %f\n", phi[i]);
          }
          // printf("xLoc : %f\n", xLoc);
          for (i = 0; i < theSpace->n; i++) {
            B[2*map[i] + 0] += jac * weight * phi[i] * f_x * xLoc;
            B[2*map[i] + 1] += jac * weight * phi[i] * f_y * xLoc;
          }
        }
      }
    }
  }
}

void femElasticityApplyDirichlet(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;

  for (int node = 0; node < theNodes->nNodes; node++) {
    femConstrainedNode *theConstrainedNode = &theProblem->constrainedNodes[node];
    if (theConstrainedNode->type == UNDEFINED)
      continue;
    femBoundaryType type = theConstrainedNode->type;

    if (type == DIRICHLET_X) {
      double value = theConstrainedNode->value1;
      femFullSystemConstrain(theSystem, 2 * node + 0, value);
    }
    if (type == DIRICHLET_Y) {
      double value = theConstrainedNode->value1;
      femFullSystemConstrain(theSystem, 2 * node + 1, value);
    }
    if (type == DIRICHLET_XY) {
      double value_x = theConstrainedNode->value1;
      double value_y = theConstrainedNode->value2;
      femFullSystemConstrain(theSystem, 2 * node + 0, value_x);
      femFullSystemConstrain(theSystem, 2 * node + 1, value_y);
    }

    if (type == DIRICHLET_N) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      femFullSystemConstrainDirichlet_N(theSystem, node, value, nx, ny);
    }
    if (type == DIRICHLET_T) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      femFullSystemConstrainDirichlet_T(theSystem, node, value, nx, ny);
    }
    if (type == DIRICHLET_NT) {
      double value_n = theConstrainedNode->value1;
      double value_t = theConstrainedNode->value2;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      femFullSystemConstrainDirichlet_N(theSystem, node, value_n, nx, ny);
      femFullSystemConstrainDirichlet_T(theSystem, node, value_t, nx, ny);
    }
  }
}

double *femElasticitySolve(femProblem *theProblem) {
  femElasticityAssembleElements(theProblem);
  femElasticityAssembleNeumann(theProblem);
  femElasticityApplyDirichlet(theProblem);
  double *soluce = solve_cg(theProblem->system);  // solveur itératif des gradients conjugués
  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}




