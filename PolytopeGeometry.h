#ifndef POLYTOPEGEOMETRY_H
#define POLYTOPEGEOMETRY_H

//#include "CommonDef.h"
#include "PyramidTree.h"
#include <array>
#include <vector>
#include <memory>


// Maximum number of components for domain and range coordinates.
const unsigned int MaxRangeDim = 10;
const unsigned int MaxDomainDim = 10;

// Maximum number of slices per polytope.
const unsigned int MaxNrSlicesPerPtope = 512;

// Bounds on polytope structure.
const int MaxNrFacetsPerPtope = 50;
const int MaxNrPointsPerPtope = 5000;


// Record structure for active polytopes.
struct _Ptope {
    IdType start;        // offset of first facet id in Facet array.
    short size;             // number of facets
    short dim;              // spatial dimensionality
    int center;             // if non-null, the center of the polytope.
    int minCoord;           // if non-null, the lowest corner of bounding box
    int splitDim;           // dimension on which any cached fragments were generated
    int fragments;    // if non-null, offset of cached fragments
    int ends;         // if non-null, offset of cached ends from cutting planes
    };


class PolytopeGeometry {
public:

  // Description:
  // Initialize a geometry object, specifying domain and range dimensions,
  // and providing in advance.
  PolytopeGeometry(const short m_topo, const short m_geom, const short n, const double *domainBounds);  
  ~PolytopeGeometry(); 
  // Description:
  // Clear the active polytope store to empty, and initialise it to
  // receive new polytopes within the specified domain (spatial)
  // bounds.  "domainBounds" points to 2M doubles containing the
  // min & max bounds for each of the M ordinates.
  void ResetForNextCell(const double *domainBounds);

  // Returns the number of stored Ptopes
  IdType  GetNrStoredPtopes(){return NrStoredPtopes; };

  // Return the number of stored facets  
  IdType  GetFacetsNr(){return CenterFacets.size(); };

  //Return a read only handle to CenterFaces 
  const IdType * GetCenterFacets(){return &CenterFacets[0][0];};
  void  ShdFacetSids(IdType id, IdType& uid, IdType& vid){uid = CenterFacets[id][0]; vid = CenterFacets[id][1];};

  void GetFacetCenter(IdType id, double * p);
  
  // Description:
  // Specialized version of PrintSelf, which prints the top
  // level (polytope id and list of points, domain and range
  // coordinates) for the "sz" active polytopes whose indices
  // are pointed to bt "top".
//?  void PrintToplevel(ostream &os, vtkIdType *top, int sz);

  // Description:
  // Copy information on active polytopes into the store.
  // Active polytopes are not changed/deleted.
   void StoreActivePolytopes(IdType *ids, int nr);

  // Description:
  // Add a simplex specified by a list of point ids, and
  // a dimension, to the active polytopes.  The point ids
  // must refer to coordinates that have already been added.
  IdType BuildFromSimplex(IdType *pt, int d);

  // Description:
  // Add a vertex (0-d polytope), specified by a coordinate in 
  // domain space, and a coordinate in range space.
  IdType AddCoord(const double *dc,const double *rc);

  // Description:
  // Add a polytope to the active polytope set.  The polytope
  // is specified by a dimensionality, and a list of "sz" ids
  // to the polytope/vertex ids of its facets.
  IdType AddPtope(IdType *pt, int d, int sz);

  // Description:
  // Access the domain coordinate of an active polytope,
  // returning coordinate in buffer provided by caller.
  void GetActiveDomCoord(IdType id, double *c);

  // Description:
  // Access the range coordinate of an active polytope,
  // returning coordinate in buffer provided by caller.
  void GetActiveRngCoord(IdType id, double *c);

  // Description:
  // Return a specific component of the domain coordinate
  // of an active polytope.
//?  double GetActiveDomCoordComponent(vtkIdType id, int c);

  // Description:
  // Return a specific component of the range coordinate
  // of an active polytope.
  double GetActiveRngCoordComponent(IdType id, int c);

  // Description:
  // Return a specific component of the range coordinate
  // of a stored polytope.
//?  double GetStoredRangeComponent(vtkIdType id, int c);
  double  GetStoredRangeComponent(IdType id, int c);

  // Description:
  // Split the active polytope with id "ptid" into a set 
  // of polytope fragments by slicing ptid against "nrts" 
  // threshold values.  Each threshold value specifies a
  // cutting plane orthogonal to axis "field"; thresholds
  // values must appear in increasing value.
  // Each intersections of a cutting plane and the polytope
  // also result in a lower-dimensional 'end' facet.  The
  // fragments and ends are added to the active polytopes, 
  // and their ids returned in the caller-allocated 
  // "frags" and "ends" buffers.  
  // *NOTE*: the buffers must
  // contain at least one more slot than the number of
  // cutting planes: the ith slot in frags holds the 
  // fragment that is LESS THAN the ith cutting plane;
  // the "nrts" slot holds the residue of the polytope 
  // that is greater than the highest cutting plane.
  void Split(
    IdType ptid, 
    double *thresholds, 
    int nrts, 
    int field, 
    IdType *frags, 
    IdType *ends);

/*  // Description:
  // Retrieve the object holding the center point of
  // each stored polytope facet.
  //vtkGetObjectMacro(CenterPoints, vtkPoints);
  vtkGetObjectMacro(CenterLocator,vtkPyramidTree);

  // Description:
  // Retrieve the two-component array giving the ids
  // of the one or two polytopes that share a facet.  
  vtkGetObjectMacro(CenterFacets, vtkIdTypeArray);

  // Description:
  // Return the number of stored polytopes.
  vtkGetMacro(NrStoredPtopes, vtkIdType);
*/  
 void PrintActivePtope(IdType id);  
private:
//?  vtkPolytopeGeometry(const vtkPolytopeGeometry&); // Not implemented
//?  void operator=(const vtkPolytopeGeometry&); // Not implemented

  static const int TABSIZE  =      9001;
  static const int SETSIZE  =      32;
  static const int MaxNrPtopes =    9999;
  static const int MaxNrFacets =  50000;
  
  // dimensionality
  short M_topo;  // domain (spatial) dimension of the topology
  short M_geom;  // domain (spatial) dimension of the geometry
  short N;       // range (data) dimension

  // Active polytope storage.  
  // Polytopes are held in an array of records indexed directly by
  // polytope ID.  For algorithms that require fast lookup of polytopes
  // based on structure, a hashtable is used, hashing on the *SET* of
  // points defining the polytope.  The hash key is computed from the
  // *sorted* list of polytope components (points or facet polytope ids).
  // For efficiency (avoiding repeated memory allocation) the facet ids
  // of all polytopes are held in a single extensible array, with each
  // polytope holding a start offset and number of facets.
  struct _Ptope Ptope[MaxNrPtopes];
  IdType Facet[MaxNrFacets];

  // Hash table:
  // We use open hashing, implemented by three structures:
  // - Index records the polytope id for a given hash.
  // - SetSize records the key size (i.e. number of polytope components)
  //   for occupied entries, and -1 for unused entries.  It thuse doubles
  //   as an occupancy map and fast reject test for matches against a key.
  // - PointSet holds the key for each hash entry.  
  //   CURRENTLY ASSUMES FIXED UPPER-BOUND (SETSIZE) ON KEY SIZE.

  IdType Index[TABSIZE];
  short SetSize[TABSIZE];
  IdType PointSet[TABSIZE][SETSIZE];

  // Global counters for active polytope and facet tables.
  IdType NextFacet;
  IdType NextPtope;

  // index for stored polytopes
  IdType NrStoredPtopes;


  // Helper functions for Split method, handling
  // specific cases.
  void SplitEdge(
    struct _Ptope *entry, 
    double *thresholds, 
    int nrts, 
    int field,
    IdType *frags, 
    IdType *ends);
    
  void SplitPtope(
    struct _Ptope *entry, 
    double *thresholds, 
    int nrts, 
    int field,
    IdType *frags, 
    IdType *ends);


  // Helper function for Split method, extracting a polytope
  // of a specific dimensionality from candidate set either by
  // finding a single instance of the required size, or by
  // assembling a candidate from a set of lower-dimensional facets.
  IdType Select(
    IdType *pts, 
    int nrts,
    int nrps,
    int slice,
    int startDim, int endDim);


  double *GetPtopeCoordMin(IdType id);
//?  void GetPtopeFacets(vtkIdType id, vtkIdType *&ids, int &nr);
  double *GetPtopeCenter(IdType id);
  
  void EnsureCacheCapacity(IdType size); 
  void EnsureDblCacheCapacity(IdType size);


  // Acceleraton strutures for active polytopes
  std::unique_ptr<PyramidTree> ActiveDomCoord;
     

  // Storage and accelerated lookup for facet center points.
  std::unique_ptr<PyramidTree> CenterLocator;


  // Map linking polytopes to facet centers.
  // Each m-dimensional polytope has a set of 
  // (m-1)-dimensional facets.  We assume each facet 
  // is shared by at most two polytopes.  The two 
  // components in the ith entry of CenterFacets are 
  // the point ids of the one (or two) polytopes 
  // containing the ith in the CenterPoints array
  // as defined above.
  std::vector< std::array<IdType, 2> > CenterFacets;
  IdType *RecCache;
  IdType MaxCache;
  IdType CacheSize;

  double *DblCache;
  IdType MaxDblCache;
  IdType DblCacheSize;


  // Storage for range coordinates (multi-field scalars).
  std::vector< std::vector<double>> ActiveScalars;
  std::vector< std::vector<double>> StoredScalars;
  
};

#endif 
