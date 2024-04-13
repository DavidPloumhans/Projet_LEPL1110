#include "fem.h"


void femElasticityAssembleElements(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];
        } 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; 
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; 
            }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                            dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                            dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                            dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                            dphidx[i] * c * dphidx[j]) * jac * weight; 
                }
            }
            // voir pour B, mais à priori c'est bon dans notre cas
        }
    }
}

void femElasticityAssembleNeumann(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->ruleEdge;
    femDiscrete    *theSpace = theProblem->spaceEdge;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[2],y[2],phi[2];
    int iBnd,iElem,iInteg,iEdge,i,j,d,map[2],mapU[2];
    int nLocal = 2;
    double *B  = theSystem->B;

    for(iBnd=0; iBnd < theProblem->nBoundaryConditions; iBnd++){
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];  // theCondition est la iBnd-ème condition
        femBoundaryType type = theCondition->type;  // Soit NEUMANN_X, NEUMANN_Y, DIRICHLET_X ou DIRICHLET_Y
        double value = theCondition->value;

        // Traite que NEUMANN_X et NEUMANN_Y ici
        // calcule l'intégrale avec la règle d'intégration
        for (iEdge = 0; iEdge < theEdges->nElem; iEdge++) {  // Pour chaque arête

            int shift=-1;
            if (type == DIRICHLET_X)  shift = 0;      
            if (type == DIRICHLET_Y)  shift = 1;  
            if (shift == -1) return; 
            int *elem = theCondition->domain->elem;
            int nElem = theCondition->domain->nElem;

            for (j=0; j < nLocal; j++) {
                // mapper les bons noeuds en faisant le lien de theEdges->elem au noeuds globaux n'est pas trivial
                for (i=0; i<2; i++) {
                    int node = theCondition->domain->mesh->elem[2*elem[j]+i];
                }
            }
            
            double jac = 0.5 * sqrt((x[1]-x[0])*(x[1]-x[0]) + (y[1]-y[0])*(y[1]-y[0]));  // Le Jacobien est la moitié de la longeur du segement 1D
            for (iInteg=0; iInteg < theRule->n; iInteg++) {  // Pour chaque point d'inégration
                double xsi    = theRule->xsi[iInteg];
                double weight = theRule->weight[iInteg];  
                femDiscretePhi(theSpace,xsi,phi);  // une seule fonction de forme en 1D donc pas besoin de phi2 ni de eta. Eta n'est pas donné par la règle d'intégration et c'est pour ça que t'avais un crash
                
                for (i = 0; i < theSpace->n; i++) {
                    B[mapU[i]] += phi[i] * value * jac * weight; 
                }
            }
        }
        
    }
}

double *femElasticitySolve(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
    femFullSystemInit(theSystem);
    femElasticityAssembleElements(theProblem);    
    femElasticityAssembleNeumann(theProblem);
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); 
        }
    }
    femFullSystemEliminate(theSystem);  // faudra le changer et mettre le solveur bande
    memcpy(theProblem->soluce, theSystem->B, theSystem->size * sizeof(double));
    return theProblem->soluce;
}