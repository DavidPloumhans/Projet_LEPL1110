/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

femGeo theGeometry;

femGeo *geoGetGeometry(void) { return &theGeometry; }

double geoSizeDefault(double x, double y) { return theGeometry.h; }

double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data) { return theGeometry.geoSize(x, y); }
void geoInitialize(void) {
  int ierr;
  theGeometry.geoSize = geoSizeDefault;
  gmshInitialize(0, NULL, 1, 0, &ierr);
  ErrorGmsh(ierr);
  gmshModelAdd("MyGeometry", &ierr);
  ErrorGmsh(ierr);
  gmshModelMeshSetSizeCallback(geoGmshSize, NULL, &ierr);
  ErrorGmsh(ierr);
  theGeometry.theNodes = NULL;
  theGeometry.theElements = NULL;
  theGeometry.theEdges = NULL;
  theGeometry.nDomains = 0;
  theGeometry.theDomains = NULL;
}

void geoFree(void) {
  if (theGeometry.theNodes) {
    free(theGeometry.theNodes->X);
    free(theGeometry.theNodes->Y);
    free(theGeometry.theNodes->number);
    free(theGeometry.theNodes);
  }
  if (theGeometry.theElements) {
    free(theGeometry.theElements->elem);
    free(theGeometry.theElements);
  }
  if (theGeometry.theEdges) {
    free(theGeometry.theEdges->elem);
    free(theGeometry.theEdges);
  }
  for (int i = 0; i < theGeometry.nDomains; i++) {
    free(theGeometry.theDomains[i]->elem);
    free(theGeometry.theDomains[i]);
  }
  free(theGeometry.theDomains);
}

void geoFinalize(void) {
  int ierr;
  geoFree();
  gmshFinalize(&ierr);
  ErrorGmsh(ierr);
}

void geoSetSizeCallback(double (*geoSize)(double x, double y)) { theGeometry.geoSize = geoSize; }

void geoMeshImport(void) {
  int ierr;

  /* Importing nodes */

  size_t nNode, n, m, *node;
  double *xyz, *trash;
  //  gmshModelMeshRenumberNodes(&ierr);                        ErrorGmsh(ierr);
  gmshModelMeshGetNodes(&node, &nNode, &xyz, &n, &trash, &m, -1, -1, 0, 0, &ierr);
  ErrorGmsh(ierr);
  femNodes *theNodes = malloc(sizeof(femNodes));
  theNodes->nNodes = nNode;
  theNodes->X = malloc(sizeof(double) * (theNodes->nNodes));
  theNodes->Y = malloc(sizeof(double) * (theNodes->nNodes));
  for (int i = 0; i < theNodes->nNodes; i++) {
    theNodes->X[i] = xyz[3 * node[i] - 3];
    theNodes->Y[i] = xyz[3 * node[i] - 2];
  }
  theGeometry.theNodes = theNodes;
  gmshFree(node);
  gmshFree(xyz);
  gmshFree(trash);
  printf("Geo     : Importing %d nodes \n", theGeometry.theNodes->nNodes);

  /* Importing elements */
  /* Pas super joli : a ameliorer pour eviter la triple copie */

  // Triangles
  size_t nElem, *elem;
  gmshModelMeshGetElementsByType(2, &elem, &nElem, &node, &nNode, -1, 0, 1, &ierr);
  ErrorGmsh(ierr);
  if (nElem != 0) {
    femMesh *theElements = malloc(sizeof(femMesh));
    theElements->nLocalNode = 3;
    theElements->nodes = theNodes;
    theElements->nElem = nElem;
    theElements->elem = malloc(sizeof(int) * 3 * theElements->nElem);
    for (int i = 0; i < theElements->nElem; i++)
      for (int j = 0; j < theElements->nLocalNode; j++)
        theElements->elem[3 * i + j] = node[3 * i + j] - 1;
    theGeometry.theElements = theElements;
    gmshFree(node);
    gmshFree(elem);
    printf("Geo     : Importing %d triangles \n", theElements->nElem);
  }

  // Quads
  int nElemTriangles = nElem;
  gmshModelMeshGetElementsByType(3, &elem, &nElem, &node, &nNode, -1, 0, 1, &ierr);
  ErrorGmsh(ierr);
  if (nElem != 0 && nElemTriangles != 0)
    Error("Cannot consider hybrid geometry with triangles and quads :-(");

  if (nElem != 0) {
    femMesh *theElements = malloc(sizeof(femMesh));
    theElements->nLocalNode = 4;
    theElements->nodes = theNodes;
    theElements->nElem = nElem;
    theElements->elem = malloc(sizeof(int) * 4 * theElements->nElem);
    for (int i = 0; i < theElements->nElem; i++)
      for (int j = 0; j < theElements->nLocalNode; j++)
        theElements->elem[4 * i + j] = node[4 * i + j] - 1;
    theGeometry.theElements = theElements;
    gmshFree(node);
    gmshFree(elem);
    printf("Geo     : Importing %d quads \n", theElements->nElem);
  }

  // Compute node renumbering
  femMesh *theElements = theGeometry.theElements;
  int *connectedNodes = calloc(theNodes->nNodes, sizeof(int));
  for (int iElem = 0; iElem < theElements->nElem; iElem++) {
    for (int i = 0; i < theElements->nLocalNode; i++) {
      connectedNodes[theElements->elem[iElem * theElements->nLocalNode + i]] = 1;
    }
  }
  int *nodeRenumber = malloc(theNodes->nNodes * sizeof(int));
  int countNodes = 0;
  for (int i = 0; i < theNodes->nNodes; i++) {
    if (connectedNodes[i]) {
      nodeRenumber[i] = countNodes;
      countNodes++;
    } else {
      nodeRenumber[i] = -2147483648;
    }
  }

  // condensing nodes
  for (int i = 0; i < theNodes->nNodes; i++) {
    if (nodeRenumber[i] < 0)
      continue;
    theNodes->X[nodeRenumber[i]] = theNodes->X[i];
    theNodes->Y[nodeRenumber[i]] = theNodes->Y[i];
  }
  theNodes->nNodes = countNodes;
  theNodes->X = realloc(theNodes->X, sizeof(double) * (theNodes->nNodes));
  theNodes->Y = realloc(theNodes->Y, sizeof(double) * (theNodes->nNodes));

  // renumbering elements
  int nLocalNode = theElements->nLocalNode;
  for (int i = 0; i < theElements->nElem; i++) {
    for (int j = 0; j < nLocalNode; j++) {
      theElements->elem[nLocalNode * i + j] = nodeRenumber[theElements->elem[nLocalNode * i + j]];
    }
  }

  gmshModelMeshGetElementsByType(1, &elem, &nElem, &node, &nNode, -1, 0, 1, &ierr);
  ErrorGmsh(ierr);
  femMesh *theEdges = malloc(sizeof(femMesh));
  theEdges->nLocalNode = 2;
  theEdges->nElem = nElem;
  theEdges->nodes = theNodes;
  theEdges->elem = malloc(sizeof(int) * 2 * theEdges->nElem);
  int countEdges = 0;
  int *connectedEdges = calloc(nElem, sizeof(int));
  int *edgeRenumber = malloc(nElem * sizeof(int));

  for (int i = 0; i < nElem; i++) {
    int map[2] = {node[2 * i + 0] - 1, node[2 * i + 1] - 1};
    if (!connectedNodes[map[0]] || !connectedNodes[map[1]]) {
      edgeRenumber[i] = -2147483648;
      continue;
    }
    for (int j = 0; j < theEdges->nLocalNode; j++) {
      connectedEdges[i] = 1;
      theEdges->elem[2 * countEdges + j] = nodeRenumber[map[j]];
    }
    edgeRenumber[i] = countEdges;
    countEdges++;
  }
  theEdges->nElem = countEdges;
  theEdges->elem = realloc(theEdges->elem, sizeof(int) * 2 * theEdges->nElem);
  theGeometry.theEdges = theEdges;
  int shiftEdges = elem[0];
  gmshFree(node);
  gmshFree(elem);
  printf("Geo     : Importing %d edges \n", theEdges->nElem);

  /* Importing 1D entities */
  int *dimTags;
  gmshModelGetEntities(&dimTags, &n, 1, &ierr);
  ErrorGmsh(ierr);
  theGeometry.nDomains = n / 2;
  theGeometry.theDomains = malloc(sizeof(femDomain *) * n / 2);
  printf("Geo     : Importing %d entities \n", theGeometry.nDomains);

  for (int i = 0; i < n / 2; i++) {
    int dim = dimTags[2 * i + 0];
    int tag = dimTags[2 * i + 1];
    femDomain *theDomain = malloc(sizeof(femDomain));
    theGeometry.theDomains[i] = theDomain;
    theDomain->mesh = theEdges;
    sprintf(theDomain->name, "Entity %d ", tag - 1);

    int *elementType;
    size_t nElementType, **elementTags, *nElementTags, nnElementTags, **nodesTags, *nNodesTags, nnNodesTags;
    gmshModelMeshGetElements(&elementType, &nElementType, &elementTags, &nElementTags, &nnElementTags, &nodesTags, &nNodesTags, &nnNodesTags, dim, tag, &ierr);
    theDomain->nElem = nElementTags[0];
    theDomain->elem = malloc(sizeof(int) * 2 * theDomain->nElem);
    int nElemCount = 0;
    for (int j = 0; j < theDomain->nElem; j++) {
      int mapEdge = edgeRenumber[elementTags[0][j] - shiftEdges];
      if (mapEdge < 0)
        continue;
      theDomain->elem[nElemCount] = mapEdge;
      nElemCount++;
    }
    theDomain->nElem = nElemCount;
    printf("Geo     : Entity %d : %d elements \n", i, theDomain->nElem);
    gmshFree(nElementTags);
    gmshFree(nNodesTags);
    gmshFree(elementTags);
    gmshFree(nodesTags);
    gmshFree(elementType);
  }
  gmshFree(dimTags);

  free(connectedNodes);
  free(nodeRenumber);
  free(connectedEdges);
  free(edgeRenumber);

  // filter out empty domains
  int countDomains = 0;
  for (int i = 0; i < theGeometry.nDomains; i++) {
    if (theGeometry.theDomains[i]->nElem != 0) {
      theGeometry.theDomains[countDomains] = theGeometry.theDomains[i];
      countDomains++;
    } else {
      free(theGeometry.theDomains[i]->elem);
      free(theGeometry.theDomains[i]);
      theGeometry.theDomains[i] = NULL;
    }
  }
  theGeometry.nDomains = countDomains;

  return;
}

void geoMeshPrint(void) {
  femNodes *theNodes = theGeometry.theNodes;
  if (theNodes != NULL) {
    printf("Number of nodes %d \n", theNodes->nNodes);
    for (int i = 0; i < theNodes->nNodes; i++) {
      printf("%6d : %14.7e %14.7e \n", i, theNodes->X[i], theNodes->Y[i]);
    }
  }
  femMesh *theEdges = theGeometry.theEdges;
  if (theEdges != NULL) {
    printf("Number of edges %d \n", theEdges->nElem);
    int *elem = theEdges->elem;
    for (int i = 0; i < theEdges->nElem; i++) {
      printf("%6d : %6d %6d \n", i, elem[2 * i], elem[2 * i + 1]);
    }
  }
  femMesh *theElements = theGeometry.theElements;
  if (theElements != NULL) {
    if (theElements->nLocalNode == 3) {
      printf("Number of triangles %d \n", theElements->nElem);
      int *elem = theElements->elem;
      for (int i = 0; i < theElements->nElem; i++) {
        printf("%6d : %6d %6d %6d\n", i, elem[3 * i], elem[3 * i + 1], elem[3 * i + 2]);
      }
    }
    if (theElements->nLocalNode == 4) {
      printf("Number of quads %d \n", theElements->nElem);
      int *elem = theElements->elem;
      for (int i = 0; i < theElements->nElem; i++) {
        printf("%6d : %6d %6d %6d %6d\n", i, elem[4 * i], elem[4 * i + 1], elem[4 * i + 2], elem[4 * i + 3]);
      }
    }
  }
  int nDomains = theGeometry.nDomains;
  printf("Number of domains %d\n", nDomains);
  for (int iDomain = 0; iDomain < nDomains; iDomain++) {
    femDomain *theDomain = theGeometry.theDomains[iDomain];
    printf("  Domain : %6d \n", iDomain);
    printf("  Name : %s\n", theDomain->name);
    printf("  Number of elements : %6d\n", theDomain->nElem);
    for (int i = 0; i < theDomain->nElem; i++) {
      //         if (i != theDomain->nElem  && (i % 10) != 0)  printf(" - ");
      printf("%6d", theDomain->elem[i]);
      if ((i + 1) != theDomain->nElem && (i + 1) % 10 == 0)
        printf("\n");
    }
    printf("\n");
  }
}

void geoMeshWrite(const char *filename) {
  FILE *file = fopen(filename, "w");
  if (!file) {
    printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename);
    exit(-1);
  }

  femNodes *theNodes = theGeometry.theNodes;
  fprintf(file, "Number of nodes %d \n", theNodes->nNodes);
  for (int i = 0; i < theNodes->nNodes; i++) {
    fprintf(file, "%6d : %14.7e %14.7e \n", i, theNodes->X[i], theNodes->Y[i]);
  }

  femMesh *theEdges = theGeometry.theEdges;
  fprintf(file, "Number of edges %d \n", theEdges->nElem);
  int *elem = theEdges->elem;
  for (int i = 0; i < theEdges->nElem; i++) {
    fprintf(file, "%6d : %6d %6d \n", i, elem[2 * i], elem[2 * i + 1]);
  }

  femMesh *theElements = theGeometry.theElements;
  if (theElements->nLocalNode == 3) {
    fprintf(file, "Number of triangles %d \n", theElements->nElem);
    elem = theElements->elem;
    for (int i = 0; i < theElements->nElem; i++) {
      fprintf(file, "%6d : %6d %6d %6d\n", i, elem[3 * i], elem[3 * i + 1], elem[3 * i + 2]);
    }
  }
  if (theElements->nLocalNode == 4) {
    fprintf(file, "Number of quads %d \n", theElements->nElem);
    elem = theElements->elem;
    for (int i = 0; i < theElements->nElem; i++) {
      fprintf(file, "%6d : %6d %6d %6d %6d\n", i, elem[4 * i], elem[4 * i + 1], elem[4 * i + 2], elem[4 * i + 3]);
    }
  }

  int nDomains = theGeometry.nDomains;
  fprintf(file, "Number of domains %d\n", nDomains);
  for (int iDomain = 0; iDomain < nDomains; iDomain++) {
    femDomain *theDomain = theGeometry.theDomains[iDomain];
    fprintf(file, "  Domain : %6d \n", iDomain);
    fprintf(file, "  Name : %s\n", theDomain->name);
    fprintf(file, "  Number of elements : %6d\n", theDomain->nElem);
    for (int i = 0; i < theDomain->nElem; i++) {
      fprintf(file, "%6d", theDomain->elem[i]);
      if ((i + 1) != theDomain->nElem && (i + 1) % 10 == 0)
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
  }

  fclose(file);
}

void geoMeshRead(const char *filename) {
  FILE *file = fopen(filename, "r");
  if (!file) {
    printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename);
    exit(-1);
  }

  int trash, *elem;

  femNodes *theNodes = malloc(sizeof(femNodes));
  theGeometry.theNodes = theNodes;
  ErrorScan(fscanf(file, "Number of nodes %d \n", &theNodes->nNodes));
  theNodes->X = malloc(sizeof(double) * (theNodes->nNodes));
  theNodes->Y = malloc(sizeof(double) * (theNodes->nNodes));
  theNodes->number = malloc(sizeof(int) * (theNodes->nNodes));
  for (int i = 0; i < theNodes->nNodes; i++) {
    ErrorScan(fscanf(file, "%d : %le %le \n", &theNodes->number[i], &theNodes->X[i], &theNodes->Y[i]));
  }

  femMesh *theEdges = malloc(sizeof(femMesh));
  theGeometry.theEdges = theEdges;
  theEdges->nLocalNode = 2;
  theEdges->nodes = theNodes;
  ErrorScan(fscanf(file, "Number of edges %d \n", &theEdges->nElem));
  theEdges->elem = malloc(sizeof(int) * theEdges->nLocalNode * theEdges->nElem);
  for (int i = 0; i < theEdges->nElem; ++i) {
    elem = theEdges->elem;
    ErrorScan(fscanf(file, "%6d : %6d %6d \n", &trash, &elem[2 * i], &elem[2 * i + 1]));
  }

  femMesh *theElements = malloc(sizeof(femMesh));
  theGeometry.theElements = theElements;
  theElements->nLocalNode = 0;
  theElements->nodes = theNodes;
  char elementType[MAXNAME];
  ErrorScan(fscanf(file, "Number of %s %d \n", elementType, &theElements->nElem));
  if (strncasecmp(elementType, "triangles", MAXNAME) == 0) {
    theElements->nLocalNode = 3;
    theElements->elem = malloc(sizeof(int) * theElements->nLocalNode * theElements->nElem);
    for (int i = 0; i < theElements->nElem; ++i) {
      elem = theElements->elem;
      ErrorScan(fscanf(file, "%6d : %6d %6d %6d \n", &trash, &elem[3 * i], &elem[3 * i + 1], &elem[3 * i + 2]));
    }
  }
  if (strncasecmp(elementType, "quads", MAXNAME) == 0) {
    theElements->nLocalNode = 4;
    theElements->elem = malloc(sizeof(int) * theElements->nLocalNode * theElements->nElem);
    for (int i = 0; i < theElements->nElem; ++i) {
      elem = theElements->elem;
      ErrorScan(fscanf(file, "%6d : %6d %6d %6d %6d \n", &trash, &elem[4 * i], &elem[4 * i + 1], &elem[4 * i + 2], &elem[4 * i + 3]));
    }
  }

  ErrorScan(fscanf(file, "Number of domains %d\n", &theGeometry.nDomains));
  int nDomains = theGeometry.nDomains;
  theGeometry.theDomains = malloc(sizeof(femDomain *) * nDomains);
  for (int iDomain = 0; iDomain < nDomains; iDomain++) {
    femDomain *theDomain = malloc(sizeof(femDomain));
    theGeometry.theDomains[iDomain] = theDomain;
    theDomain->mesh = theEdges;
    ErrorScan(fscanf(file, "  Domain : %6d \n", &trash));
    ErrorScan(fscanf(file, "  Name : %[^\n]s \n", (char *)&theDomain->name));
    ErrorScan(fscanf(file, "  Number of elements : %6d\n", &theDomain->nElem));
    theDomain->elem = malloc(sizeof(int) * 2 * theDomain->nElem);
    for (int i = 0; i < theDomain->nElem; i++) {
      ErrorScan(fscanf(file, "%6d", &theDomain->elem[i]));
      if ((i + 1) != theDomain->nElem && (i + 1) % 10 == 0)
        ErrorScan(fscanf(file, "\n"));
    }
  }

  fclose(file);
}

void geoSetDomainName(int iDomain, char *name) {
  if (iDomain >= theGeometry.nDomains)
    Error("Illegal domain number");
  if (geoGetDomain(name) != -1)
    Error("Cannot use the same name for two domains");
  sprintf(theGeometry.theDomains[iDomain]->name, "%s", name);
}

int geoGetDomain(char *name) {
  int theIndex = -1;
  int nDomains = theGeometry.nDomains;
  for (int iDomain = 0; iDomain < nDomains; iDomain++) {
    femDomain *theDomain = theGeometry.theDomains[iDomain];
    if (strncasecmp(name, theDomain->name, MAXNAME) == 0)
      theIndex = iDomain;
  }
  return theIndex;
}

static const double _gaussQuad4Xsi[4] = {-0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Eta[4] = {0.577350269189626, -0.577350269189626, -0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Weight[4] = {1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000};
static const double _gaussTri3Xsi[3] = {0.166666666666667, 0.666666666666667, 0.166666666666667};
static const double _gaussTri3Eta[3] = {0.166666666666667, 0.166666666666667, 0.666666666666667};
static const double _gaussTri3Weight[3] = {0.166666666666667, 0.166666666666667, 0.166666666666667};
static const double _gaussEdge2Xsi[2] = {0.577350269189626, -0.577350269189626};
static const double _gaussEdge2Weight[2] = {1.000000000000000, 1.000000000000000};

femIntegration *femIntegrationCreate(int n, femElementType type) {
  femIntegration *theRule = malloc(sizeof(femIntegration));
  if (type == FEM_QUAD && n == 4) {
    theRule->n = 4;
    theRule->xsi = _gaussQuad4Xsi;
    theRule->eta = _gaussQuad4Eta;
    theRule->weight = _gaussQuad4Weight;
  } else if (type == FEM_TRIANGLE && n == 3) {
    theRule->n = 3;
    theRule->xsi = _gaussTri3Xsi;
    theRule->eta = _gaussTri3Eta;
    theRule->weight = _gaussTri3Weight;
  } else if (type == FEM_EDGE && n == 2) {
    theRule->n = 2;
    theRule->xsi = _gaussEdge2Xsi;
    theRule->eta = NULL;
    theRule->weight = _gaussEdge2Weight;
  } else
    Error("Cannot create such an integration rule !");
  return theRule;
}

void femIntegrationFree(femIntegration *theRule) { free(theRule); }

void _q1c0_x(double *xsi, double *eta) {
  xsi[0] = 1.0;
  eta[0] = 1.0;
  xsi[1] = -1.0;
  eta[1] = 1.0;
  xsi[2] = -1.0;
  eta[2] = -1.0;
  xsi[3] = 1.0;
  eta[3] = -1.0;
}

void _q1c0_phi(double xsi, double eta, double *phi) {
  phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;
  phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
  phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
  phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
}

void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta) {
  dphidxsi[0] = (1.0 + eta) / 4.0;
  dphidxsi[1] = -(1.0 + eta) / 4.0;
  dphidxsi[2] = -(1.0 - eta) / 4.0;
  dphidxsi[3] = (1.0 - eta) / 4.0;
  dphideta[0] = (1.0 + xsi) / 4.0;
  dphideta[1] = (1.0 - xsi) / 4.0;
  dphideta[2] = -(1.0 - xsi) / 4.0;
  dphideta[3] = -(1.0 + xsi) / 4.0;
}

void _p1c0_x(double *xsi, double *eta) {
  xsi[0] = 0.0;
  eta[0] = 0.0;
  xsi[1] = 1.0;
  eta[1] = 0.0;
  xsi[2] = 0.0;
  eta[2] = 1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi) {
  phi[0] = 1 - xsi - eta;
  phi[1] = xsi;
  phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta) {
  dphidxsi[0] = -1.0;
  dphidxsi[1] = 1.0;
  dphidxsi[2] = 0.0;
  dphideta[0] = -1.0;
  dphideta[1] = 0.0;
  dphideta[2] = 1.0;
}

void _e1c0_x(double *xsi) {
  xsi[0] = -1.0;
  xsi[1] = 1.0;
}

void _e1c0_phi(double xsi, double *phi) {
  phi[0] = (1 - xsi) / 2.0;
  phi[1] = (1 + xsi) / 2.0;
}

void _e1c0_dphidx(double xsi, double *dphidxsi) {
  dphidxsi[0] = -0.5;
  dphidxsi[1] = 0.5;
}

femDiscrete *femDiscreteCreate(int n, femElementType type) {
  femDiscrete *theSpace = malloc(sizeof(femDiscrete));
  if (type == FEM_QUAD && n == 4) {
    theSpace->n = 4;
    theSpace->x2 = _q1c0_x;
    theSpace->phi2 = _q1c0_phi;
    theSpace->dphi2dx = _q1c0_dphidx;
  } else if (type == FEM_TRIANGLE && n == 3) {
    theSpace->n = 3;
    theSpace->x2 = _p1c0_x;
    theSpace->phi2 = _p1c0_phi;
    theSpace->dphi2dx = _p1c0_dphidx;
  } else if (type == FEM_EDGE && n == 2) {
    theSpace->n = 2;
    theSpace->x = _e1c0_x;
    theSpace->phi = _e1c0_phi;
    theSpace->dphidx = _e1c0_dphidx;
  } else
    Error("Cannot create such a discrete space !");
  return theSpace;
}

void femDiscreteFree(femDiscrete *theSpace) { free(theSpace); }

void femDiscreteXsi(femDiscrete *mySpace, double *xsi) { mySpace->x(xsi); }

void femDiscretePhi(femDiscrete *mySpace, double xsi, double *phi) { mySpace->phi(xsi, phi); }

void femDiscreteXsi2(femDiscrete *mySpace, double *xsi, double *eta) { mySpace->x2(xsi, eta); }

void femDiscretePhi2(femDiscrete *mySpace, double xsi, double eta, double *phi) { mySpace->phi2(xsi, eta, phi); }

void femDiscreteDphi2(femDiscrete *mySpace, double xsi, double eta, double *dphidxsi, double *dphideta) { mySpace->dphi2dx(xsi, eta, dphidxsi, dphideta); }

void femDiscretePrint(femDiscrete *mySpace) {
  int i, j;
  int n = mySpace->n;
  double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];

  femDiscreteXsi2(mySpace, xsi, eta);
  for (i = 0; i < n; i++) {

    femDiscretePhi2(mySpace, xsi[i], eta[i], phi);
    femDiscreteDphi2(mySpace, xsi[i], eta[i], dphidxsi, dphideta);

    for (j = 0; j < n; j++) {
      printf("(xsi=%+.1f,eta=%+.1f) : ", xsi[i], eta[i]);
      printf(" phi(%d)=%+.1f", j, phi[j]);
      printf("   dphidxsi(%d)=%+.1f", j, dphidxsi[j]);
      printf("   dphideta(%d)=%+.1f \n", j, dphideta[j]);
    }
    printf(" \n");
  }
}

femFullSystem *femFullSystemCreate(int size) {
  femFullSystem *theSystem = malloc(sizeof(femFullSystem));
  femFullSystemAlloc(theSystem, size);
  femFullSystemInit(theSystem);

  return theSystem;
}

void femFullSystemFree(femFullSystem *theSystem) {
  free(theSystem->A);
  free(theSystem->B);
  free(theSystem);
}

void femFullSystemAlloc(femFullSystem *mySystem, int size) {
  int i;
  double *elem = malloc(sizeof(double) * size * (size + 1));
  mySystem->A = malloc(sizeof(double *) * size);
  mySystem->B = elem;
  mySystem->A[0] = elem + size;
  mySystem->size = size;
  for (i = 1; i < size; i++)
    mySystem->A[i] = mySystem->A[i - 1] + size;
}

void femFullSystemInit(femFullSystem *mySystem) {
  int i, size = mySystem->size;
  for (i = 0; i < size * (size + 1); i++)
    mySystem->B[i] = 0;
}

void femFullSystemPrint(femFullSystem *mySystem) {
  double **A, *B;
  int i, j, size;

  A = mySystem->A;
  B = mySystem->B;
  size = mySystem->size;

  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++)
      if (A[i][j] == 0)
        printf("         ");
      else
        printf(" %+.1e", A[i][j]);
    printf(" :  %+.1e \n", B[i]);
  }
}

double *femFullSystemEliminate(femFullSystem *mySystem) {
  double **A, *B, factor;
  int i, j, k, size;

  A = mySystem->A;
  B = mySystem->B;
  size = mySystem->size;

  /* Gauss elimination */

  for (k = 0; k < size; k++) {
    if (fabs(A[k][k]) <= 1e-16) {
      printf("Pivot index %d  ", k);
      printf("Pivot value %e  ", A[k][k]);
      Error("Cannot eliminate with such a pivot");
    }
    for (i = k + 1; i < size; i++) {
      factor = A[i][k] / A[k][k];
      for (j = k + 1; j < size; j++)
        A[i][j] = A[i][j] - A[k][j] * factor;
      B[i] = B[i] - B[k] * factor;
    }
  }

  /* Back-substitution */

  for (i = size - 1; i >= 0; i--) {
    factor = 0;
    for (j = i + 1; j < size; j++)
      factor += A[i][j] * B[j];
    B[i] = (B[i] - factor) / A[i][i];
  }

  return (mySystem->B);
}

void femFullSystemConstrain(femFullSystem *mySystem, int myNode, double myValue) {
  double **A, *B;
  int i, size;

  A = mySystem->A;
  B = mySystem->B;
  size = mySystem->size;

  for (i = 0; i < size; i++) {
    B[i] -= myValue * A[i][myNode];
    A[i][myNode] = 0;
  }

  for (i = 0; i < size; i++)
    A[myNode][i] = 0;

  A[myNode][myNode] = 1;
  B[myNode] = myValue;
}

femProblem *femElasticityCreate(femGeo *theGeometry, double E, double nu, double rho, double gx, double gy, femElasticCase iCase) {
  femProblem *theProblem = malloc(sizeof(femProblem));
  theProblem->E = E;
  theProblem->nu = nu;
  theProblem->gx = gx;
  theProblem->gy = gy;
  theProblem->rho = rho;

  if (iCase == PLANAR_STRESS) {
    theProblem->A = E / (1 - nu * nu);
    theProblem->B = E * nu / (1 - nu * nu);
    theProblem->C = E / (2 * (1 + nu));
  } else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
    theProblem->A = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
    theProblem->B = E * nu / ((1 + nu) * (1 - 2 * nu));
    theProblem->C = E / (2 * (1 + nu));
  }

  theProblem->planarStrainStress = iCase;
  theProblem->nBoundaryConditions = 0;
  theProblem->conditions = NULL;

  int nNodes = theGeometry->theNodes->nNodes;
  int size = 2 * nNodes;
  theProblem->constrainedNodes = malloc(nNodes * sizeof(femConstrainedNode));
  for (int i = 0; i < nNodes; i++) {
    theProblem->constrainedNodes[i].type = UNDEFINED;
    theProblem->constrainedNodes[i].nx = NAN;
    theProblem->constrainedNodes[i].ny = NAN;
    theProblem->constrainedNodes[i].value2 = NAN;
    theProblem->constrainedNodes[i].value2 = NAN;
  }

  theProblem->geometry = theGeometry;
  if (theGeometry->theElements->nLocalNode == 3) {
    theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
    theProblem->rule = femIntegrationCreate(3, FEM_TRIANGLE);
  }
  if (theGeometry->theElements->nLocalNode == 4) {
    theProblem->space = femDiscreteCreate(4, FEM_QUAD);
    theProblem->rule = femIntegrationCreate(4, FEM_QUAD);
  }
  theProblem->spaceEdge = femDiscreteCreate(2, FEM_EDGE);
  theProblem->ruleEdge = femIntegrationCreate(2, FEM_EDGE);
  theProblem->system = femFullSystemCreate(size);
  return theProblem;
}

void femElasticityFree(femProblem *theProblem) {
  femFullSystemFree(theProblem->system);
  femIntegrationFree(theProblem->rule);
  femDiscreteFree(theProblem->space);
  for (int i = 0; i < theProblem->nBoundaryConditions; i++)
    free(theProblem->conditions[i]);
  free(theProblem->conditions);
  free(theProblem->soluce);
  free(theProblem->residuals);
  free(theProblem->constrainedNodes);
  free(theProblem);
}

/*
`value2` is only used for `DIRICHLET_XY` and `DIRICHLET_NT` boundary conditions. Otherwise it is ignored and set to NAN.
*/
void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value1, double value2) {
  int iDomain = geoGetDomain(nameDomain);
  if (iDomain == -1)
    Error("Undefined domain :-(");
  value2 = ((type != DIRICHLET_XY) && (type != DIRICHLET_NT)) ? NAN : value2;

  femBoundaryCondition *theBoundary = malloc(sizeof(femBoundaryCondition));
  theBoundary->domain = theProblem->geometry->theDomains[iDomain];
  theBoundary->value1 = value1;
  theBoundary->value2 = value2;
  theBoundary->type = type;
  theProblem->nBoundaryConditions++;
  int nBoundaryConditions = theProblem->nBoundaryConditions;

  if (theProblem->conditions == NULL) {
    theProblem->conditions = malloc(nBoundaryConditions * sizeof(femBoundaryCondition *));
  }

  femNodes *theNodes = theProblem->geometry->theNodes;
  if (type == DIRICHLET_X || type == DIRICHLET_Y || type == DIRICHLET_XY || type == DIRICHLET_N || type == DIRICHLET_T || type == DIRICHLET_NT) {
    // Ensure that there is only one Dirichlet boundary condition per domain
    for (int i = 0; i < nBoundaryConditions - 1; i++) {
      if (theProblem->conditions[i]->domain != theBoundary->domain)
        continue;
      femBoundaryType type_i = theProblem->conditions[i]->type;
      if (type_i == DIRICHLET_X || type_i == DIRICHLET_Y || type_i == DIRICHLET_XY || type_i == DIRICHLET_N || type_i == DIRICHLET_T || type_i == DIRICHLET_NT) {
        printf("\nTrying to set a second Dirichlet boundary condition on domain \"%s\"", nameDomain);
        Error("Only one Dirichlet boundary condition is allowed per domain");
      }
    }

    femDomain *theDomain = theProblem->geometry->theDomains[iDomain];
    int *elem = theDomain->elem;
    int nElem = theDomain->nElem;
    femConstrainedNode constrainedNode;
    constrainedNode.type = type;
    constrainedNode.value1 = value1;
    constrainedNode.value2 = value2;
    constrainedNode.nx = NAN;
    constrainedNode.ny = NAN;
    if (type == DIRICHLET_X || type == DIRICHLET_Y || type == DIRICHLET_XY) {
      for (int iElem = 0; iElem < nElem; iElem++) {
        for (int i = 0; i < 2; i++) {
          int node = theDomain->mesh->elem[2 * elem[iElem] + i];
          theProblem->constrainedNodes[node] = constrainedNode;
        }
      }
    } else { // need to compute normals
      int nNodes = theNodes->nNodes;
      double *NX = malloc(nNodes * sizeof(double));
      double *NY = malloc(nNodes * sizeof(double));
      for (int iElem = 0; iElem < nElem; iElem++) {
        int node0 = theDomain->mesh->elem[2 * elem[iElem] + 0];
        int node1 = theDomain->mesh->elem[2 * elem[iElem] + 1];
        NX[node0] = 0;
        NY[node0] = 0;
        NX[node1] = 0;
        NY[node1] = 0;
      }

      for (int iElem = 0; iElem < nElem; iElem++) {
        int node0 = theDomain->mesh->elem[2 * elem[iElem] + 0];
        int node1 = theDomain->mesh->elem[2 * elem[iElem] + 1];
        double tx = theNodes->X[node1] - theNodes->X[node0];
        double ty = theNodes->Y[node1] - theNodes->Y[node0];
        double nx = ty;
        double ny = -tx;
        NX[node0] += nx;
        NY[node0] += ny;
        NX[node1] += nx;
        NY[node1] += ny;
      }

      for (int iElem = 0; iElem < nElem; iElem++) {
        for (int i = 0; i < 2; i++) {
          int node = theDomain->mesh->elem[2 * elem[iElem] + i];
          double nx = NX[node];
          double ny = NY[node];
          double norm = hypot(nx, ny);
          theProblem->constrainedNodes[node] = constrainedNode;
          theProblem->constrainedNodes[node].nx = nx / norm;
          theProblem->constrainedNodes[node].ny = ny / norm;
        }
      }
      free(NX);
      free(NY);
    }
  }

  theProblem->conditions = realloc(theProblem->conditions, nBoundaryConditions * sizeof(femBoundaryCondition *));
  theProblem->conditions[nBoundaryConditions - 1] = theBoundary;
}

void femElasticityPrint(femProblem *theProblem) {
  printf("\n\n ======================================================================================= \n\n");
  printf(" Linear elasticity problem \n");
  printf("   Young modulus   E   = %14.7e [N/m2]\n", theProblem->E);
  printf("   Poisson's ratio nu  = %14.7e [-]\n", theProblem->nu);
  printf("   Density         rho = %14.7e [kg/m3]\n", theProblem->rho);
  printf("   Gravity-X       gx  = %14.7e [m/s2]\n", theProblem->gx);
  printf("   Gravity-Y       gy  = %14.7e [m/s2]\n", theProblem->gy);

  if (theProblem->planarStrainStress == PLANAR_STRAIN)
    printf("   Planar strains formulation \n");
  if (theProblem->planarStrainStress == PLANAR_STRESS)
    printf("   Planar stresses formulation \n");
  if (theProblem->planarStrainStress == AXISYM)
    printf("   Axisymmetric formulation \n");

  printf("   Boundary conditions : \n");
  for (int i = 0; i < theProblem->nBoundaryConditions; i++) {
    femBoundaryCondition *theCondition = theProblem->conditions[i];
    double value1 = theCondition->value1;
    double value2 = theCondition->value2;
    printf("  %20s :", theCondition->domain->name);
    if (theCondition->type == DIRICHLET_X)
      printf(" imposing %9.2e as the horizontal displacement  \n", value1);
    if (theCondition->type == DIRICHLET_Y)
      printf(" imposing %9.2e as the vertical displacement  \n", value1);
    if (theCondition->type == DIRICHLET_XY)
      printf(" imposing %9.2e, %9.2e as the displacement  \n", value1, value2);
    if (theCondition->type == DIRICHLET_N)
      printf(" imposing %9.2e as the normal displacement  \n", value1);
    if (theCondition->type == DIRICHLET_T)
      printf(" imposing %9.2e as the tangential displacement  \n", value1);
    if (theCondition->type == DIRICHLET_NT)
      printf(" imposing (%9.2e, %9.2e) as the normal and tangential displacement  \n", value1, value2);
    if (theCondition->type == NEUMANN_X)
      printf(" imposing %9.2e as the horizontal force  \n", value1);
    if (theCondition->type == NEUMANN_Y)
      printf(" imposing %9.2e as the vertical force  \n", value1);
    if (theCondition->type == NEUMANN_N)
      printf(" imposing %9.2e as the normal force  \n", value1);
    if (theCondition->type == NEUMANN_T)
      printf(" imposing %9.2e as the tangential force  \n", value1);
  }
  printf(" ======================================================================================= \n\n");
}

void femElasticityWrite(femProblem *theProblem, const char *filename) {
  FILE *file = fopen(filename, "w");
  if (!file) {
    printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename);
    exit(-1);
  }

  switch (theProblem->planarStrainStress) {
  case PLANAR_STRESS:
    fprintf(file, "Type of problem    :  Planar stresses  \n");
    break;
  case PLANAR_STRAIN:
    fprintf(file, "Type of problem    :  Planar strains \n");
    break;
  case AXISYM:
    fprintf(file, "Type of problem    :  Axi-symetric problem \n");
    break;
  default:
    fprintf(file, "Type of problem    :  Undefined  \n");
    break;
  }
  fprintf(file, "Young modulus      : %14.7e  \n", theProblem->E);
  fprintf(file, "Poisson ratio      : %14.7e  \n", theProblem->nu);
  fprintf(file, "Mass density       : %14.7e  \n", theProblem->rho);
  fprintf(file, "Gravity-X          : %14.7e  \n", theProblem->gx);
  fprintf(file, "Gravity-Y          : %14.7e  \n", theProblem->gy);

  for (int i = 0; i < theProblem->nBoundaryConditions; i++) {
    femBoundaryCondition *theCondition = theProblem->conditions[i];
    double value1 = theCondition->value1;
    double value2 = theCondition->value2;
    fprintf(file, "Boundary condition : ");
    switch (theCondition->type) {
    case DIRICHLET_X:
      fprintf(file, " Dirichlet-X        = %14.7e, %14.7e ", value1, NAN);
      break;
    case DIRICHLET_Y:
      fprintf(file, " Dirichlet-Y        = %14.7e, %14.7e ", value1, NAN);
      break;
    case DIRICHLET_XY:
      fprintf(file, " Dirichlet-XY       = %14.7e, %14.7e ", value1, value2);
      break;
    case DIRICHLET_N:
      fprintf(file, " Dirichlet-N        = %14.7e, %14.7e ", value1, NAN);
      break;
    case DIRICHLET_T:
      fprintf(file, " Dirichlet-T        = %14.7e, %14.7e ", value1, NAN);
      break;
    case DIRICHLET_NT:
      fprintf(file, " Dirichlet-NT       = %14.7e, %14.7e ", value1, value2);
      break;
    case NEUMANN_X:
      fprintf(file, " Neumann-X          = %14.7e, %14.7e ", value1, NAN);
      break;
    case NEUMANN_Y:
      fprintf(file, " Neumann-Y          = %14.7e, %14.7e ", value1, NAN);
      break;
    case NEUMANN_N:
      fprintf(file, " Neumann-N          = %14.7e, %14.7e ", value1, NAN);
      break;
    case NEUMANN_T:
      fprintf(file, " Neumann-T          = %14.7e, %14.7e ", value1, NAN);
      break;
    default:
      fprintf(file, " Undefined          = %14.7e, %14.7e ", NAN, NAN);
      break;
    }

    fprintf(file, ": %s\n", theCondition->domain->name);
  }
  fclose(file);
}

femProblem *femElasticityRead(femGeo *theGeometry, const char *filename) {
  FILE *file = fopen(filename, "r");
  if (!file) {
    printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename);
    exit(-1);
  }
  femProblem *theProblem = malloc(sizeof(femProblem));
  theProblem->nBoundaryConditions = 0;
  theProblem->conditions = NULL;

  int nNodes = theGeometry->theNodes->nNodes;
  int size = 2 * nNodes;
  theProblem->soluce = malloc(size * sizeof(double));
  theProblem->residuals = malloc(size * sizeof(double));
  for (int i = 0; i < size; i++) {
    theProblem->soluce[i] = 0.0;
    theProblem->residuals[i] = 0.0;
  }

  theProblem->constrainedNodes = malloc(nNodes * sizeof(femConstrainedNode));
  for (int i = 0; i < nNodes; i++) {
    theProblem->constrainedNodes[i].type = UNDEFINED;
    theProblem->constrainedNodes[i].nx = NAN;
    theProblem->constrainedNodes[i].ny = NAN;
    theProblem->constrainedNodes[i].value2 = NAN;
    theProblem->constrainedNodes[i].value2 = NAN;
  }

  theProblem->geometry = theGeometry;
  if (theGeometry->theElements->nLocalNode == 3) {
    theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
    theProblem->rule = femIntegrationCreate(3, FEM_TRIANGLE);
  }
  if (theGeometry->theElements->nLocalNode == 4) {
    theProblem->space = femDiscreteCreate(4, FEM_QUAD);
    theProblem->rule = femIntegrationCreate(4, FEM_QUAD);
  }
  theProblem->spaceEdge = femDiscreteCreate(2, FEM_EDGE);
  theProblem->ruleEdge = femIntegrationCreate(2, FEM_EDGE);
  theProblem->system = femFullSystemCreate(nNodes);  // petit test

  char theLine[MAXNAME];
  char theDomain[MAXNAME];
  char theArgument[MAXNAME];
  double value1, value2;
  femBoundaryType typeCondition;

  while (!feof(file)) {
    ErrorScan(fscanf(file, "%19[^\n]s \n", (char *)&theLine));
    if (strncasecmp(theLine, "Type of problem     ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %[^\n]s \n", (char *)&theArgument));
      if (strncasecmp(theArgument, "Planar stresses", 13) == 0)
        theProblem->planarStrainStress = PLANAR_STRESS;
      if (strncasecmp(theArgument, "Planar strains", 13) == 0)
        theProblem->planarStrainStress = PLANAR_STRAIN;
      if (strncasecmp(theArgument, "Axi-symetric problem", 13) == 0)
        theProblem->planarStrainStress = AXISYM;
    }
    if (strncasecmp(theLine, "Young modulus       ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %le\n", &theProblem->E));
    }
    if (strncasecmp(theLine, "Poisson ratio       ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %le\n", &theProblem->nu));
    }
    if (strncasecmp(theLine, "Mass density        ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %le\n", &theProblem->rho));
    }
    if (strncasecmp(theLine, "Gravity-X           ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %le\n", &theProblem->gx));
    }
    if (strncasecmp(theLine, "Gravity-Y           ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %le\n", &theProblem->gy));
    }
    if (strncasecmp(theLine, "Boundary condition  ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %19s = %le, %le : %[^\n]s\n", (char *)&theArgument, &value1, &value2, (char *)&theDomain));
      if (strncasecmp(theArgument, "Dirichlet-X", 19) == 0)
        typeCondition = DIRICHLET_X;
      if (strncasecmp(theArgument, "Dirichlet-Y", 19) == 0)
        typeCondition = DIRICHLET_Y;
      if (strncasecmp(theArgument, "Dirichlet-XY", 19) == 0)
        typeCondition = DIRICHLET_XY;
      if (strncasecmp(theArgument, "Dirichlet-N", 19) == 0)
        typeCondition = DIRICHLET_N;
      if (strncasecmp(theArgument, "Dirichlet-T", 19) == 0)
        typeCondition = DIRICHLET_T;
      if (strncasecmp(theArgument, "Dirichlet-NT", 19) == 0)
        typeCondition = DIRICHLET_NT;
      if (strncasecmp(theArgument, "Neumann-X", 19) == 0)
        typeCondition = NEUMANN_X;
      if (strncasecmp(theArgument, "Neumann-Y", 19) == 0)
        typeCondition = NEUMANN_Y;
      if (strncasecmp(theArgument, "Neumann-N", 19) == 0)
        typeCondition = NEUMANN_N;
      if (strncasecmp(theArgument, "Neumann-T", 19) == 0)
        typeCondition = NEUMANN_T;
      femElasticityAddBoundaryCondition(theProblem, theDomain, typeCondition, value1, value2);
    }
    ErrorScan(fscanf(file, "\n"));
  }

  int iCase = theProblem->planarStrainStress;
  double E = theProblem->E;
  double nu = theProblem->nu;

  if (iCase == PLANAR_STRESS) {
    theProblem->A = E / (1 - nu * nu);
    theProblem->B = E * nu / (1 - nu * nu);
    theProblem->C = E / (2 * (1 + nu));
  } else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
    theProblem->A = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
    theProblem->B = E * nu / ((1 + nu) * (1 - 2 * nu));
    theProblem->C = E / (2 * (1 + nu));
  }

  fclose(file);
  return theProblem;
}

femProblem *femElasticityRead2(femGeo *theGeometry, const char *filename) {
  FILE *file = fopen(filename, "r");
  if (!file) {
    printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename);
    exit(-1);
  }
  femProblem *theProblem = malloc(sizeof(femProblem));
  theProblem->nBoundaryConditions = 0;
  theProblem->conditions = NULL;

  int nNodes = theGeometry->theNodes->nNodes;
  int size = 2 * nNodes;
  theProblem->soluce = malloc(size * sizeof(double));
  theProblem->residuals = malloc(size * sizeof(double));
  for (int i = 0; i < size; i++) {
    theProblem->soluce[i] = 0.0;
    theProblem->residuals[i] = 0.0;
  }

  theProblem->constrainedNodes = malloc(nNodes * sizeof(femConstrainedNode));
  for (int i = 0; i < nNodes; i++) {
    theProblem->constrainedNodes[i].type = UNDEFINED;
    theProblem->constrainedNodes[i].nx = NAN;
    theProblem->constrainedNodes[i].ny = NAN;
    theProblem->constrainedNodes[i].value2 = NAN;
    theProblem->constrainedNodes[i].value2 = NAN;
  }

  theProblem->geometry = theGeometry;
  if (theGeometry->theElements->nLocalNode == 3) {
    theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
    theProblem->rule = femIntegrationCreate(3, FEM_TRIANGLE);
  }
  if (theGeometry->theElements->nLocalNode == 4) {
    theProblem->space = femDiscreteCreate(4, FEM_QUAD);
    theProblem->rule = femIntegrationCreate(4, FEM_QUAD);
  }
  theProblem->spaceEdge = femDiscreteCreate(2, FEM_EDGE);
  theProblem->ruleEdge = femIntegrationCreate(2, FEM_EDGE);
  theProblem->system = femFullSystemCreate(size);

  char theLine[MAXNAME];
  char theDomain[MAXNAME];
  char theArgument[MAXNAME];
  double value1, value2;
  femBoundaryType typeCondition;

  while (!feof(file)) {
    ErrorScan(fscanf(file, "%19[^\n]s \n", (char *)&theLine));
    if (strncasecmp(theLine, "Type of problem     ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %[^\n]s \n", (char *)&theArgument));
      if (strncasecmp(theArgument, "Planar stresses", 13) == 0)
        theProblem->planarStrainStress = PLANAR_STRESS;
      if (strncasecmp(theArgument, "Planar strains", 13) == 0)
        theProblem->planarStrainStress = PLANAR_STRAIN;
      if (strncasecmp(theArgument, "Axi-symetric problem", 13) == 0)
        theProblem->planarStrainStress = AXISYM;
    }
    if (strncasecmp(theLine, "Young modulus       ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %le\n", &theProblem->E));
    }
    if (strncasecmp(theLine, "Poisson ratio       ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %le\n", &theProblem->nu));
    }
    if (strncasecmp(theLine, "Mass density        ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %le\n", &theProblem->rho));
    }
    if (strncasecmp(theLine, "Gravity-X           ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %le\n", &theProblem->gx));
    }
    if (strncasecmp(theLine, "Gravity-Y           ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %le\n", &theProblem->gy));
    }
    if (strncasecmp(theLine, "Boundary condition  ", 19) == 0) {
      ErrorScan(fscanf(file, ":  %19s = %le, %le : %[^\n]s\n", (char *)&theArgument, &value1, &value2, (char *)&theDomain));
      if (strncasecmp(theArgument, "Dirichlet-X", 19) == 0)
        typeCondition = DIRICHLET_X;
      if (strncasecmp(theArgument, "Dirichlet-Y", 19) == 0)
        typeCondition = DIRICHLET_Y;
      if (strncasecmp(theArgument, "Dirichlet-XY", 19) == 0)
        typeCondition = DIRICHLET_XY;
      if (strncasecmp(theArgument, "Dirichlet-N", 19) == 0)
        typeCondition = DIRICHLET_N;
      if (strncasecmp(theArgument, "Dirichlet-T", 19) == 0)
        typeCondition = DIRICHLET_T;
      if (strncasecmp(theArgument, "Dirichlet-NT", 19) == 0)
        typeCondition = DIRICHLET_NT;
      if (strncasecmp(theArgument, "Neumann-X", 19) == 0)
        typeCondition = NEUMANN_X;
      if (strncasecmp(theArgument, "Neumann-Y", 19) == 0)
        typeCondition = NEUMANN_Y;
      if (strncasecmp(theArgument, "Neumann-N", 19) == 0)
        typeCondition = NEUMANN_N;
      if (strncasecmp(theArgument, "Neumann-T", 19) == 0)
        typeCondition = NEUMANN_T;
      femElasticityAddBoundaryCondition(theProblem, theDomain, typeCondition, value1, value2);
    }
    ErrorScan(fscanf(file, "\n"));
  }

  int iCase = theProblem->planarStrainStress;
  double E = theProblem->E;
  double nu = theProblem->nu;

  if (iCase == PLANAR_STRESS) {
    theProblem->A = E / (1 - nu * nu);
    theProblem->B = E * nu / (1 - nu * nu);
    theProblem->C = E / (2 * (1 + nu));
  } else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
    theProblem->A = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
    theProblem->B = E * nu / ((1 + nu) * (1 - 2 * nu));
    theProblem->C = E / (2 * (1 + nu));
  }

  fclose(file);
  return theProblem;
}

void femSolutionWrite(int nNodes, int nfields, double *data, const char *filename) {
  FILE *file = fopen(filename, "w");
  if (!file) {
    printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename);
    exit(-1);
  }
  fprintf(file, "Size %d,%d\n", nNodes, nfields);
  for (int i = 0; i < nNodes; i++) {
    for (int j = 0; j < nfields - 1; j++) {
      fprintf(file, "%.18le,", data[i * nfields + j]);
    }
    fprintf(file, "%.18le", data[i * nfields + nfields - 1]);
    fprintf(file, "\n");
  }
  fclose(file);
}

int femSolutiondRead(int allocated_size, double *value, const char *filename) {
  FILE *file = fopen(filename, "r");
  if (!file) {
    printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename);
    exit(-1);
  }
  int nNodes, nFields;
  ErrorScan(fscanf(file, "Size %d,%d\n", &nNodes, &nFields));
  if (nNodes * nFields > allocated_size) {
    printf("Error: allocated size is %d, but the solution file has %d nodes and %d fields", allocated_size, nNodes, nFields);
    Error("The allocated size is too small for femSolutiondRead");
  }
  for (int i = 0; i < nNodes; i++) {
    for (int j = 0; j < nFields; j++)
      ErrorScan(fscanf(file, "%le,", &value[i * nFields + j]));
    ErrorScan(fscanf(file, "\n"));
  }
  printf("Reading solution of shape (%d,%d)\n", nNodes, nFields);
  fclose(file);
  return nNodes * nFields;
}

double femMin(double *x, int n) {
  double myMin = x[0];
  int i;
  for (i = 1; i < n; i++)
    myMin = fmin(myMin, x[i]);
  return myMin;
}

double femMax(double *x, int n) {
  double myMax = x[0];
  int i;
  for (i = 1; i < n; i++)
    myMax = fmax(myMax, x[i]);
  return myMax;
}

void femError(char *text, int line, char *file) {
  printf("\n-------------------------------------------------------------------------------- ");
  printf("\n  Error in %s:%d at line %d : \n  %s\n", file, line, line, text);
  printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
  exit(69);
}

void femErrorGmsh(int ierr, int line, char *file) {
  if (ierr == 0)
    return;
  printf("\n-------------------------------------------------------------------------------- ");
  printf("\n  Error in %s:%d at line %d : \n  error code returned by gmsh %d\n", file, line, line, ierr);
  printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
  gmshFinalize(NULL);
  exit(69);
}

void femErrorScan(int test, int line, char *file) {
  if (test >= 0)
    return;

  printf("\n-------------------------------------------------------------------------------- ");
  printf("\n  Error in fscanf or fgets in %s:%d at line %d : \n", file, line, line);
  printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
  exit(69);
}

void femWarning(char *text, int line, char *file) {
  printf("\n-------------------------------------------------------------------------------- ");
  printf("\n  Warning in %s:%d at line %d : \n  %s\n", file, line, line, text);
  printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
}

void femElasticityEpsilon(femProblem *theProblem, double *E, int wtv) {  // wtv = 0 pour E_xx, 1 pour E_yy, 2 pour E_xy et 3 pour E_ThetaTheta
  femFullSystem *theSystem = theProblem->system;
  femFullSystemInit(theSystem);
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
  printf("System size : %d\n", theSystem->size);
  for (iElem = 0; iElem < theMesh->nElem; iElem++) {
    for (j = 0; j < nLocal; j++) {
      map[j] = theMesh->elem[iElem * nLocal + j];
      mapX[j] = 2 * map[j];
      mapY[j] = 2 * map[j] + 1;
      x[j] = theNodes->X[map[j]];
      y[j] = theNodes->Y[map[j]];
    }
   
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
      for (i = 0; i < theSpace->n; i++) {
        dxdxsi += x[i] * dphidxsi[i];
        dxdeta += x[i] * dphideta[i];
        dydxsi += y[i] * dphidxsi[i];
        dydeta += y[i] * dphideta[i];
      }
      double jac = dxdxsi * dydeta - dxdeta * dydxsi;
      if (jac < 0.0 && iInteg == 0) {
        printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n");
        count++;
      }
      jac = fabs(jac);
      double dudx = 0;
      double dvdy = 0;
      double xLoc = 0;
      double dudy = 0;
      double dvdx = 0;
      double uLoc = 0;
      for (i = 0; i < theSpace->n; i++) {
        dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
        dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
        dudx += dphidx[i] * theProblem->soluce[2*map[i]];
        dvdy += dphidy[i] * theProblem->soluce[2*map[i]+1];
        dudy += dphidy[i] * theProblem->soluce[2*map[i]];
        dvdx += dphidx[i] * theProblem->soluce[2*map[i]+1];
        xLoc += x[i]*phi[i];
        uLoc += theProblem->soluce[2*map[i]]*phi[i];  // fallait pas oublier le phi[i]
      }
      for (i = 0; i < theSpace->n; i++) {
        for (j = 0; j < theSpace->n; j++) {
          A[map[i]][map[j]] += phi[i]*phi[j]*weight*jac * xLoc;
        }
      }
      for (i = 0; i < theSpace->n; i++) {
        if (wtv == 0)
          B[map[i]] += phi[i] * dudx * jac * weight * xLoc;
        else if (wtv == 1)
          B[map[i]] += phi[i] * dvdy * jac * weight * xLoc;
        else if (wtv == 2)
          B[map[i]] += (0.5) * phi[i] * (dudy + dvdx) * jac * weight * xLoc;
        else if (wtv == 3)
          B[map[i]] += phi[i] * uLoc * jac * weight;
        else 
          Error("Invalid value for wtv in femElasticityEpsilon");
      }
    }
  }
  double *solution = femFullSystemEliminate(theSystem);  // il semblerait qu'il y ait eu un problème avec le solver des gradients conjugués
  memcpy(E, solution, sizeof(double) * theNodes->nNodes);
  printf("Number of elements with negative jacobian: %d\n", count);
}


void sparseMatrixFree(sparseMatrix *sp) {
  free(sp->col);
  free(sp->rptr);
  free(sp->val);
  free(sp);
}

sparseMatrix* to_sparse(double **A, int size) {
  int nnz = 0;
  for (int i = 0; i < size; i++) {
    for(int j = 0; j < size; j++) {
      nnz += (A[i][j] != 0);
    }
  }
  int* col = malloc(nnz * sizeof(int));
  int* rptr = malloc((size + 1) * sizeof(int));
  double* val = malloc(nnz * sizeof(double));
  
  nnz = 0;
  for(int i = 0; i < size; i++) {
    rptr[i] = nnz;
    for(int j = 0; j < size; j++) {
      if(A[i][j] != 0) {
        col[nnz] = j;
        val[nnz] = A[i][j];
        nnz++;
      }
    }
  }
  rptr[size] = nnz;

  sparseMatrix* sp = malloc(sizeof(sparseMatrix));
  sp->size = size;
  sp->nnz = nnz;
  sp->col = col;
  sp->rptr = rptr;
  sp->val = val;
  return sp;
}

static inline void spmv(const sparseMatrix* sp, const double* x, double* y) {
  for(int i = 0; i < sp->size; i++) {
    double s = 0;
    for(int j = sp->rptr[i]; j < sp->rptr[i+1]; j++) {
      s += sp->val[j] * x[sp->col[j]];
    }
    y[i] = s;
  }
}

static inline void residual(const sparseMatrix* sp, const double* x, const double* b, double* r) {
  spmv(sp, x, r);
  for(int i = 0; i < sp->size; i++) {
    r[i] = b[i] - r[i];
  }
}

static inline double dot(const double* x, const double* y, int size) {
  double s = 0;
  for(int i = 0; i < size; i++) {
    s += x[i] * y[i];
  }
  return s;
}

static inline void axpy(double* x, const double* y, double a, int size) {
  for(int i = 0; i < size; i++) {
    x[i] += a * y[i];
  }
}

double *solve_cg(femFullSystem *mySystem) {
  int size = mySystem->size;
  double **A = mySystem->A;
  double *B = mySystem->B;
  //clock_t start, stop;

  //start = clock();
  sparseMatrix* sp = to_sparse(A, size);
  //stop = clock();
  //printf("Time to create sparse mat: %f ms\n", 1000 * (double)(stop - start) / CLOCKS_PER_SEC);
  //start = clock();

  int niter = 0;
  double *x = malloc(size * sizeof(double));
  double *r = malloc(size * sizeof(double));
  double *p = malloc(size * sizeof(double));
  double *Ap = malloc(size * sizeof(double));

  // Initialize x, r, and p
  for (int i = 0; i < size; i++) {
    x[i] = 0.0;
    r[i] = B[i];
    p[i] = r[i];
  }

  double alpha, beta, rr, rrNew;
  rr = dot(r, r, size);

  while (niter < size) {
    spmv(sp, p, Ap);
    alpha = rr / dot(p, Ap, size);
    axpy(x, p, alpha, size);
    axpy(r, Ap, -alpha, size);
    rrNew = dot(r, r, size);
    beta = rrNew / rr;
    for (int i = 0; i < size; i++) {
      p[i] = r[i] + beta * p[i];
    }
    rr = rrNew;
    niter++;
  }

  free(r);
  free(p);
  free(Ap);

  return x;
}

double max(double a,double b,double c){
    double res;
    if(fabs(a) > fabs(b)){
        res = a;
    }
    if(fabs(c) > res){
        res = c;
    }
    return res;
}


double lev(double a[3][3],double y[3][1],double x[3][1]){

    y[0][0] = (a[0][0] * x[0][0]) + (a[0][1] * x[1][0]) + (a[0][2] * x[2][0]);
    y[1][0] = (a[1][0] * x[0][0]) + (a[1][1] * x[1][0]) + (a[1][2] * x[2][0]);
    y[2][0] = (a[2][0] * x[0][0]) + (a[2][1] * x[1][0]) + (a[2][2] * x[2][0]);
    
    
    double l1 = max(y[0][0],y[1][0],y[2][0]);

    for(int i = 0;i<3;i++){
        x[i][0] = y[i][0] / l1;
    }

    int maxIter = 10000;
    int count = 0;
    while(count < maxIter){
        y[0][0] = (a[0][0] * x[0][0]) + (a[0][1] * x[1][0]) + (a[0][2] * x[2][0]);
        y[1][0] = (a[1][0] * x[0][0]) + (a[1][1] * x[1][0]) + (a[1][2] * x[2][0]);
        y[2][0] = (a[2][0] * x[0][0]) + (a[2][1] * x[1][0]) + (a[2][2] * x[2][0]);

        double l2 = max(y[0][0],y[1][0],y[2][0]);

        for(int i = 0;i<3;i++){
            x[i][0] = y[i][0] / l2;
        }

        double error = fabs(l2 - l1) / fabs(l2);

        l1 = l2;

        if(error < 0.0001){
            return l1;
        }
        count++;
    }
    return 0;  // pas ouf mais devrait résoudre le problème
}


void calculateEigenValues(double a[3][3], double eigenvalues[3]) {
  // src : https://github.com/b-abhinay07/ICS_CP_Attagallu/blob/main/B22AI012_B22EE004_B23EE1049_B23CS1010.c
  double y[3][1];
  double x[3][1] = {{1.0},{1.0},{1.0}};

  double largesteigenvalue = lev(a,y,x);

  double sum = (a[0][0] + a[1][1] + a[2][2]) - largesteigenvalue;
  double product = ((a[0][0] * a[1][1] * a[2][2]) - (a[0][0] * a[2][1] * 
  a[1][2]) - (a[0][1] * a[1][0] * a[2][2]) + (a[0][1] * a[2][0] * a[1][2]) + 
  (a[0][2] * a[1][0] * a[2][1]) - (a[0][2] * a[2][0] * a[1][1]))/largesteigenvalue;
  double dis = (sum * sum) - (4 * product);
  double l2 = (sum + sqrt(dis)) / 2;
  double l3 = (sum - sqrt(dis)) / 2;
  if (largesteigenvalue > l2) {
    eigenvalues[0] = largesteigenvalue;
    eigenvalues[1] = l2;
    eigenvalues[2] = l3;
  } else if (largesteigenvalue < l3) {
    eigenvalues[0] = l2;
    eigenvalues[1] = l3;
    eigenvalues[2] = largesteigenvalue;
  } else {
    eigenvalues[0] = l2;
    eigenvalues[1] = largesteigenvalue;
    eigenvalues[2] = l3;
  }
}


// Performs the Jacobi eigenvalue algorithm until the largest absolute
// off-diagonal element is smaller than Epsilon.
// The matrices S and U^T of size NxN are stored in row major memory layout.
// The runtime of the algorithm is in O(N^2) per iteration.
int Jacobi(double *S, double *UT, int N, double Epsilon) {
    // Set U^T to identity matrix
    for(int K = 0; K < N*N; ++K) {
        UT[K] = 0.0;
    }
    for(int K = 0; K < N; ++K) {
        UT[K*N + K] = 1.0;
    }
    
    for(int It = 0, I, J;; ++It) {
        // Seek largest (absolute) off-diagonal element
        // Initialize with lower right most off-diagonal
        I = N - 2;
        J = N - 1;
        for(int R = 0; R < N - 2; ++R) {
            for(int C = R + 1; C < N; ++C) {
                if(fabs(S[R*N + C]) > fabs(S[I*N + J])) {
                    I = R;
                    J = C;
                }
            }
        }
        
        // Done if largest element is smaller than epsilon
        if(fabs(S[I*N + J]) < Epsilon) {
            return It;
        }
        
        // Compute sine and cosine (pick larger root because why not)
        double Sii = S[I*N + I];
        double Sjj = S[J*N + J];
        double Sij = S[I*N + J];
        double Cot2A = (Sjj - Sii)/(2.0*Sij);
        double TanA = -Cot2A + sqrt(Cot2A*Cot2A + 1.0);
        double CosA = 1.0/sqrt(1.0 + TanA*TanA);
        double SinA = TanA*CosA;
        double CSq = CosA*CosA;
        double SSq = SinA*SinA;
        double SC = SinA*CosA;
        
        // S <- R*S*R^T
        S[I*N + I] = CSq*Sii - 2*SC*Sij + SSq*Sjj;
        S[J*N + J] = SSq*Sii + 2*SC*Sij + CSq*Sjj;
        S[I*N + J] = 0.0;
        S[J*N + I] = 0.0;
        for(int K = 0; K < N; ++K) {
            if(K != I && K != J) {
                double Sik = S[I*N + K];
                double Sjk = S[J*N + K];
                S[I*N + K] = S[K*N + I] = CosA*Sik - SinA*Sjk;
                S[J*N + K] = S[K*N + J] = SinA*Sik + CosA*Sjk;
            }
        }
        
        // U^T <- U^T*R^T (U^T will contain the eigenvectors)
        for(int K = 0; K < N; ++K) {
            double UTki = UT[K*N + I];
            double UTkj = UT[K*N + J];
            UT[K*N + I] = CosA*UTki - SinA*UTkj;
            UT[K*N + J] = CosA*UTkj + SinA*UTki;
        }
    }
}


void PrintMatrix(double *M, int N) {
    for(int I = 0; I < N; ++I) {
        for(int J = 0; J < N; ++J) {
            printf("%.9f ", M[I*N + J]);
        }
        printf("\n");
    }
    printf("\n");
}

double *femElasticityForces(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
    femFullSystemInit(theSystem);
    printf("A");
    femElasticityAssembleElements(theProblem);
    printf("B");
    femElasticityAssembleNeumann(theProblem);
    printf("C");
    femElasticityApplyDirichlet(theProblem);      
    printf("In function") ;
    double *theResidual = theProblem->residuals;
    for(int i=0; i < theSystem->size; i++) {
        theResidual[i] = 0.0;}
    for(int i=0; i < theSystem->size; i++){
        for(int j=0; j < theSystem->size; j++){
            theResidual[i] += theSystem->A[i][j] * theProblem->soluce[j]; 
        }
        theResidual[i] -= theSystem->B[i]; 
    }
    return theResidual;
}

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
  printf("AA\n");
  for (iElem = 0; iElem < theMesh->nElem; iElem++) {
    printf("AAB\n");
    for (j = 0; j < nLocal; j++) {
      map[j] = theMesh->elem[iElem * nLocal + j];
      mapX[j] = 2 * map[j];
      mapY[j] = 2 * map[j] + 1;
      x[j] = theNodes->X[map[j]];
      y[j] = theNodes->Y[map[j]];
    }
    printf("AB\n");
    
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
      printf("AC\n");
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
        xLoc += x[i]*phi[i];  // x local
      }
      double jac = dxdxsi * dydeta - dxdeta * dydxsi;
      if (jac < 0.0 && iInteg == 0) {
        printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n");
        count++;
      }
      jac = fabs(jac);
      printf("jac = %f\n", jac);
      for (i = 0; i < theSpace->n; i++) {  // Calcul des dérivées des fonctions de forme
          dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
          dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
      }
      
      printf("AD\n");
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
      printf("AE\n");
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
          int xLoc = 0.0;
          double xsi = theRule->xsi[iInteg];
          double weight = theRule->weight[iInteg];
          femDiscretePhi(theSpace, xsi, phi);
          for (i = 0; i < theSpace->n; i++) {
            xLoc += x[i]*phi[i];  // x local
          }
          for (i = 0; i < theSpace->n; i++) {
            B[2*map[i] + 0] += jac * weight * phi[i] * f_x * xLoc;  // j'en suis pas méga sur
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
void femFullSystemConstrainDirichlet_N(femFullSystem *mySystem, int myNode, double myValue, double nx, double ny) {
  double **A, *B;
  int i, size;
  //MyValue correspond à la valeur de la contrainte

  A = mySystem->A;
  B = mySystem->B;
  size = mySystem->size;
  
  double tx = -ny;
  double ty = nx;
  double b_t = tx*B[2*myNode] + ty*B[2*myNode+1];
  double a_tn = nx*(tx*A[2*myNode][2*myNode]+ty*A[2*myNode+1][2*myNode]) + ny*(tx*A[2*myNode][2*myNode+1]+ty*A[2*myNode+1][2*myNode+1]);
  double a_tt = tx*(tx*A[2*myNode][2*myNode]+ty*A[2*myNode+1][2*myNode]) + ty*(tx*A[2*myNode][2*myNode+1]+ty*A[2*myNode+1][2*myNode+1]);
  //Modification des lignes et des colonnes
  for(int i = 0; i < 2*myNode; i++) {
    A[i][2*myNode] = tx*(A[i][2*myNode]*tx + A[i][2*myNode+1]*ty);
    A[i][2*myNode+1] = ty*(A[i][2*myNode]*tx + A[i][2*myNode+1]*ty);
    B[i] -= myValue*A[i][2*myNode];
  }
  for (i=2*myNode+2; i < size; i++) {
    A[i][2*myNode] = tx*(A[i][2*myNode]*tx + A[i][2*myNode+1]*ty);
    A[i][2*myNode+1] = ty*(A[i][2*myNode]*tx + A[i][2*myNode+1]*ty);
    B[i] -= myValue*A[i][2*myNode];
  }
  for (i = 0; i < 2*myNode; i++) {
      // Prendre les deux lignes à modifier
      A[2*myNode][i] = tx*(A[2*myNode][i]*tx + A[2*myNode+1][i]*ty);
      A[2*myNode+1][i] = ty*(A[2*myNode][i]*tx + A[2*myNode+1][i]*ty);
    }
  for(i = 2*myNode+2; i < size; i++) {
    A[2*myNode][i] = tx*(A[2*myNode][i]*tx + A[2*myNode+1][i]*ty);
    A[2*myNode+1][i] = ty*(A[2*myNode][i]*tx + A[2*myNode+1][i]*ty);
  }
  // Nous allons maintenant effectuer les opérations sur la petite matrice 
  // 2x2
  A[2*myNode][2*myNode] = nx*nx + tx*tx*a_tt;
  A[2*myNode][2*myNode+1] = nx*ny + tx*ty*a_tt;
  A[2*myNode+1][2*myNode] = nx*ny + ty*tx*a_tt;
  A[2*myNode+1][2*myNode+1] = ny*ny + ty*ty*a_tt;
  B[2*myNode] = myValue*nx + tx*(b_t - myValue*a_tn);
  B[2*myNode+1] = myValue*ny + ty*(b_t - myValue*a_tn); 
}

void femFullSystemConstrainDirichlet_T(femFullSystem *mySystem, int myNode, double myValue, double nx, double ny) {

  double **A, *B;
  int i, size;
  //Myvalue correspond à la valeur de la contrainte
  double tx = -ny;
  double ty = nx;
  double b_n = nx*B[2*myNode] + ny*B[2*myNode+1];
  double a_tn = nx*(tx*A[2*myNode][2*myNode]+ty*A[2*myNode+1][2*myNode]) + ny*(tx*A[2*myNode][2*myNode+1]+ty*A[2*myNode+1][2*myNode+1]);
  double a_nn = nx*(nx*A[2*myNode][2*myNode]+ny*A[2*myNode+1][2*myNode]) + ny*(nx*A[2*myNode][2*myNode+1]+ny*A[2*myNode+1][2*myNode+1]);
  //Modification des lignes et des colonnes
  for(int i = 0; i < 2*myNode; i++) {
    A[i][2*myNode] = nx*(A[i][2*myNode]*nx + A[i][2*myNode+1]*ny);
    A[i][2*myNode+1] = ny*(A[i][2*myNode]*nx + A[i][2*myNode+1]*ny);
    B[i] -= myValue*A[i][2*myNode+1];
  }
  for (i=2*myNode+2; i < size; i++) {
    A[i][2*myNode] = tx*(A[i][2*myNode]*tx + A[i][2*myNode+1]*ty);
    A[i][2*myNode+1] = ty*(A[i][2*myNode]*tx + A[i][2*myNode+1]*ty);
    B[i] -= myValue*A[i][2*myNode+1];
  }
  for (i = 0; i < 2*myNode; i++) {
      // Prendre les deux lignes à modifier
      A[2*myNode][i] = nx*(A[2*myNode][i]*nx + A[2*myNode+1][i]*ny);
      A[2*myNode+1][i] = ny*(A[2*myNode][i]*nx + A[2*myNode+1][i]*ny);
    }
  for(i = 2*myNode+2; i < size; i++) {
    A[2*myNode][i] = nx*(A[2*myNode][i]*nx + A[2*myNode+1][i]*ny);
    A[2*myNode+1][i] = ny*(A[2*myNode][i]*nx + A[2*myNode+1][i]*ny);
  }
  // Nous allons maintenant effectuer les opérations sur la petite matrice 
  // 2x2
  A[2*myNode][2*myNode] = a_nn*nx*nx + tx*tx;
  A[2*myNode][2*myNode+1] = a_nn*nx*ny + tx*ty;
  A[2*myNode+1][2*myNode] = a_nn*nx*ny + ty*tx;
  A[2*myNode+1][2*myNode+1] = a_nn*ny*ny + ty*ty;
  B[2*myNode] = myValue*tx + nx*(b_n - myValue*a_tn);
  B[2*myNode+1] = myValue*ty + ny*(b_n - myValue*a_tn); 
}

double *femElasticitySolve(femProblem *theProblem) {
  femElasticityAssembleElements(theProblem);
  femElasticityAssembleNeumann(theProblem);
  femElasticityApplyDirichlet(theProblem);

  double *soluce = femFullSystemEliminate(theProblem->system);  // solveur itératif des gradients conjugués
  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}


void calculateAnalytic(femProblem *theProblem, double *uv, int nNodes, double p) {  // p la contrainte en haut
  double E = theProblem->E;
  double nu = theProblem->nu;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  double *X = theNodes->X;
  double *Y = theNodes->Y;
  for (int i = 0; i < nNodes; i++) {
    uv[2*i] = nu / E * p * X[i];
    uv[2*i+1] = -p / E * Y[i];
  }  
}