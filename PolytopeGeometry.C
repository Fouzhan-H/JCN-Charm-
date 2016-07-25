#include "PolytopeGeometry.h"

#include <cstring>
#include <stdio.h>

#include <algorithm>
#include <math.h>

// Tolerance used when comparing components of points.
#define EPSILON 0.0000001

// ------------------------------------------------------------------

// Helper function:
// Compare two double values subject to an error tolerance.
// Return -1, 0, 1 respectively if t is less than,
// approximately equal to, or greater than, base.
int Compare(double t, double base)
{
  if (fabs(t - base) < EPSILON)
    return 0;
  if (t < base)
    return -1;
  return 1;
}

// ------------------------------------------------------------------

PolytopeGeometry::PolytopeGeometry(
  const short m_topo, 
  const short m_geom, 
  const short n,
  const double *domainBounds)
  // Create data arrays to hold range values for both active 
  // and stored polytopes.
  : ActiveScalars (0, std::vector<double>(n)) 
  , StoredScalars (0, std::vector<double>(n))
  , CenterFacets (0)
{
  if (m_topo > m_geom) {
    WarnMsg("Polytope topology is of higher dimension than geometry.");
  }

  // Initialize maximum dimensions of polytopes & data to be stored.
  this->M_topo = m_topo;
  this->M_geom = m_geom;
  this->N = n;


  // Acceleration structures for active polytopes
  this->ActiveDomCoord  = std::unique_ptr<PyramidTree> (new PyramidTree); 
  this->ActiveDomCoord->Initialize(this->M_geom, domainBounds);
    
  
  // Initialize active polytope table.
  // Ptope entry 0 is reserved as a "null" (-1-dimensional) ptope.
  this->Ptope[0].dim   = -1;
  this->Ptope[0].start =  0;
  this->Ptope[0].size  =  0;
  this->Ptope[0].minCoord = -1;
  this->Ptope[0].fragments = -1;
  this->Ptope[0].ends = -1;
  this->NextPtope = 1;

  // Counter for stored (inactive) polytopes.  
  this->NrStoredPtopes = 0;


  // Storage for facet center, and lookup.
  this->CenterLocator = std::unique_ptr<PyramidTree> (new PyramidTree); 
  this->CenterLocator->Initialize(this->M_geom, domainBounds);


//  ActiveScalars = std::unique_ptr<std::vector<std::vector<double > > > (new std::vector<std::vector<double>>(n, std::vector<double>(n) ));
  

  this->MaxCache = 0;
  this->RecCache = NULL;
  this->CacheSize = 0;

  this->MaxDblCache = 0;
  this->DblCache = NULL;
  this->DblCacheSize = 0;
  

/*  for(int i = 1; i < MaxNrPtopes; i++)
    {
    this->Ptope[i].center = -1;
    this->Ptope[i].minCoord = -1;
    this->Ptope[i].fragments = -1;
    this->Ptope[i].ends = -1;
    }
  */
}


PolytopeGeometry::~PolytopeGeometry() 
{
  delete [] this->RecCache;
  delete [] this->DblCache;

//  cout << "RecCache size/use: " << this->MaxCache << "/" << this->CacheSize << endl;
//  cout << "DblCache size/use: " << this->MaxDblCache << "/" << this->DblCacheSize << endl;
}

// ------------------------------------------------------------------
// Clear any active polytopes, resetting active storage structures,
// and re-initialize the structures for new polytopes within the
// specified spatial (domain) bounds.
void PolytopeGeometry::ResetForNextCell(const double *domainBounds)
{

/*  struct _Ptope *entry = this->Ptope+1;
  
  // Clear polytope table.
  // Entry 0, the Nil polytope, is unaffected.
  for (int i = 1; i < this->NextPtope; i++, entry++)
    {
    // All held in global caches, so offsets not pointers.
    entry->center = -1;
    entry->minCoord = -1;
    entry->fragments = -1;
    entry->ends = -1;
    }
*/
    
  // Reset polytope and facet tables counters.
  this->NextFacet = 0;
  this->NextPtope = 1;
  
  // Reset hash-table support.
  for (int i = 1; i < PolytopeGeometry::TABSIZE; i++)
    {
    this->SetSize[i] = -1;
    }

  this -> ActiveScalars.clear();
  
  // Reset domain coordinate storage and lookup acceleration.
  this->ActiveDomCoord->Initialize(this->M_geom,domainBounds);
  this->CacheSize = 0;
  this->EnsureCacheCapacity(5000);
  this->DblCacheSize = 0;
  this->EnsureDblCacheCapacity(5000);
}

// ------------------------------------------------------------------
void PolytopeGeometry::GetFacetCenter(IdType id, double* p){
  CenterLocator-> GetPoint(id,p); 
}


// ------------------------------------------------------------------


// Copy polytope information from active to stored representation.
// "ids" contains the polytope id for the top-level of the polytopes
// to be copied, "nrPtopes" is the number of top-level ids to be copied.

void PolytopeGeometry::StoreActivePolytopes(IdType *ids, int nrPtopes)
{
  IdType pid;
  bool unseen;
  struct _Ptope *entry;
  IdType facet;
  
  double *p;
  double div;
  double sc, minF, sw;
   
  // Process each top-level polytope in turn.
  for (int i = 0; i < nrPtopes; i++)
    {
    // entry refers to new top-level ptope.    
    entry = this->Ptope + ids[i];    
    p = this->GetPtopeCoordMin(ids[i]);
   
    // Copy polytope min coordinate into stored table.
    this->StoredScalars.push_back(std::vector<double> (N)); 
    this->StoredScalars[this -> NrStoredPtopes + i].assign(p, p+N); 

    // For each facet of the polytope, obtain the facet center and insert into
    // table of stored centers.  Update CenterFacets table to link the polytope
    // to the facet center id; at most two polytopes share a facet, if this was
    // the first occurrence of the facet center point, initialize the second
    // polytope slot to -1, otherwise use the second slot
        //stops at i = 1 and j = 0
    for (int j = 0; j < entry->size; j++)
      {

      facet = this->Facet[entry->start + j];
      p = this->GetPtopeCenter(facet);

      
      unseen = this->CenterLocator->InsertUniquePoint(p, pid);
          
      if (unseen)
        {
        // First occurrence of this center, use component 0 
        // and initialize the second component to -1.
        this->CenterFacets.push_back({ this->NrStoredPtopes+ i, -1});
        }
      else
        {
        // Second occurrence, use component 1.
        this->CenterFacets[pid][1] = this->NrStoredPtopes + i;
        }
      }
      
    }
  this->NrStoredPtopes += nrPtopes;

} // StoreActivePolytopes.


// ------------------------------------------------------------------

// Construct an active polytope from a set of coordinates defining
// a d-dimensional simplex.

IdType PolytopeGeometry::BuildFromSimplex(IdType *pt, int d)
{
  IdType points[MaxNrFacetsPerPtope];
  IdType facets[MaxNrFacetsPerPtope];
  IdType lines[3];
  IdType triangles[4];
  
  int size = d+1;
  IdType sid;

  switch(d) {
    
    case 0:     // must be a coordinate, already in table.
      break;
    
    case 1:     // edge defined by two points, use directly.
      sid = this->AddPtope(pt, d, d+1);
      break;

    case 2:     // triangle: build three edge polytopes, then define a polygon.
      points[0] = pt[0]; points[1] = pt[1]; lines[0] = this->AddPtope(points, 1, 2);
      points[0] = pt[0]; points[1] = pt[2]; lines[1] = this->AddPtope(points, 1, 2);
      points[0] = pt[1]; points[1] = pt[2]; lines[2] = this->AddPtope(points, 1, 2);
      sid = this->AddPtope(lines, 2, 3);
      break;
      
    case 3:     // tetrahedron: build four triangles, one per face.
      // Face 0 1 2
      points[0] = pt[0]; points[1] = pt[1]; lines[0] = this->AddPtope(points, 1, 2);
      points[0] = pt[0]; points[1] = pt[2]; lines[1] = this->AddPtope(points, 1, 2);
      points[0] = pt[1]; points[1] = pt[2]; lines[2] = this->AddPtope(points, 1, 2);
      triangles[0] = this->AddPtope(lines, 2, 3);
    
      // Face 0 1 3
      points[0] = pt[0]; points[1] = pt[1]; lines[0] = this->AddPtope(points, 1, 2);
      points[0] = pt[0]; points[1] = pt[3]; lines[1] = this->AddPtope(points, 1, 2);
      points[0] = pt[1]; points[1] = pt[3]; lines[2] = this->AddPtope(points, 1, 2);
      triangles[1] = this->AddPtope(lines, 2, 3);

      // Face 0 2 3
      points[0] = pt[0]; points[1] = pt[2]; lines[0] = this->AddPtope(points, 1, 2);
      points[0] = pt[0]; points[1] = pt[3]; lines[1] = this->AddPtope(points, 1, 2);
      points[0] = pt[2]; points[1] = pt[3]; lines[2] = this->AddPtope(points, 1, 2);
      triangles[2] = this->AddPtope(lines, 2, 3);

      // Face 1 2 3
      points[0] = pt[1]; points[1] = pt[2]; lines[0] = this->AddPtope(points, 1, 2);
      points[0] = pt[1]; points[1] = pt[3]; lines[1] = this->AddPtope(points, 1, 2);
      points[0] = pt[2]; points[1] = pt[3]; lines[2] = this->AddPtope(points, 1, 2);
      triangles[3] = this->AddPtope(lines, 2, 3);
      sid = this->AddPtope(triangles, 3, 4);
      break;
      
    default:
      // Reduce dimensionality by one, and recurse.
      // Given n points, form n (d-1)-dimensional facets by
      // systematically omitting one point at a time.
      for (int omit = 0; omit < size; omit++)
        {
        for (int i = 0, p = 0; i < size; i++)
          {
          if (i == omit) continue;
          points[p++] = pt[i];
          }
        facets[omit] = this->BuildFromSimplex(points, d-1);
        } // for each face
      sid = this->AddPtope(facets, d, d+1);        

     } // switch dimension

  // return the id of the top-level polytope.
  return sid;
}

// ---------------------------------------------------------------------

// Add a vertex on a polytope; the vertex is defined by 
// a domain coordinate (dc) and a range coordinate (rc).

IdType PolytopeGeometry::AddCoord(const double *dc, const double *rc)
{ 
  IdType pid;
  int unseen;
  double *base;

  // The *domain* coordinate is unique to each polytope vertex;
  // check if the coordinate is already known, if so return the
  // existing coordinate id.
  unseen = this->ActiveDomCoord->InsertUniquePoint(dc, pid);
  // if this is a *new* domain point, store the corresponding range.
  if (unseen){
    this->ActiveScalars.push_back( std::vector<double> (this->N) );
    this->ActiveScalars[pid].assign(rc, rc + this->N);  
  }

  // Coordinate ids are represented as negative numbers.
  // Encode pid as a coordinate and return.
  return -1 - pid;
}

// ---------------------------------------------------------------------

// Add a polytope specified by a list of 'sz' point/sub-tope indices.

IdType PolytopeGeometry::AddPtope(IdType *pt, int d, int sz)
{
  IdType p;
  IdType temp[MaxNrPointsPerPtope];
  unsigned long long h;
  struct _Ptope *entry;
  double fp, fq;
  double c1[MaxRangeDim], c2[MaxRangeDim];

  IdType a, b, c; 

  // Copy component ids to avoid interference with caller's data, and
  // sort ids to provide a unique hash key for a given set of ids.
  for (int i = 0; i < sz; i++) 
    {
    temp[i] = pt[i];
    }
  std::sort(temp, temp+sz);

  // Hash the key.  Hash function to use is determined by number of indices
  // in the polytope, to try and get a wider spread of values over the common
  // low-dimensional cases.
   switch (sz) {
    case 1: 
      h = temp[0]; 
      break;
    case 2:
      h = (temp[1] << 16) + temp[0];
      break;
    case 3:
      h = (temp[2] << 24) + (temp[1] << 16) + temp[0];
      break;
    default:
      h = 0;
      for (int j = 0; j < sz; j++) { h += temp[j]; h = h << 4; }
    }

  // compute index into hash table.
  
  h = h % PolytopeGeometry::TABSIZE;

  // linear probing: starting from initial hash location, probe
  // consecutive locations until either an empty slot is found,
  // or an exact match is found.  In the former case, create 
  // and initialize a new polytope.

  while (true)
    {
    if (this->SetSize[h] < 0)
      {
      // Empty slot, poltope not present.
      // Copy polytope key into this slot.
      this->SetSize[h] = sz;
      for (int i = 0; i < sz; i++)
        {
        this->PointSet[h][i] = temp[i];
        }
      // Initialize the polytope record.
      entry = this->Ptope + this->NextPtope;
      entry->start = this->NextFacet;
      entry->dim   = d;
      entry->size  = sz;
      entry->center = -1;
      entry->fragments = -1;
      entry->ends = -1;
      // Copy facet indices into facet array,
      for (int i = 0; i < sz; i++)
        {
        this->Facet[this->NextFacet++] = temp[i];
        }
      
       this->EnsureDblCacheCapacity(this->DblCacheSize + this->N);
       entry->minCoord = this->DblCacheSize;
       this->DblCacheSize += this->N;
      
      // Compute minimum coordinate of polytope.
      if (d == 1)
        {
        // Edge: minimum is componentwise min over 
        // coordinates of end points.
        this->GetActiveRngCoord(pt[0], c1);
        this->GetActiveRngCoord(pt[1], c2);
        for (int j = 0; j < this->N; j++)
          {
          this->DblCache[entry->minCoord + j] = std::min(c1[j], c2[j]);
          }
        }
      else
        {
        // Non-edge: Initialise minimum to min coordinate of
        // first component, then update against min of each
        // subsequent component.
        for (int j = 0; j < this->N; j++)
          {
          this->DblCache[entry->minCoord + j] = this->DblCache[this->Ptope[pt[0]].minCoord + j];
          }
        for (int i = 1; i < sz; i++)
          {
          for (int j = 0; j < this->N; j++)
            {
            this->DblCache[entry->minCoord + j] = std::min( 
                this->DblCache[entry->minCoord + j],
                this->DblCache[this->Ptope[pt[i]].minCoord + j] );
            }
          }
        }
      
      // Update hash index entry to point to this polytope.
      this->Index[h] = this->NextPtope++;
      break;
      } // empty slot 
    else 
      {
      // Hash table contains an entry at this position.
      // Quick test: does the size match?
      if (this->SetSize[h] == sz)
        {
        // Test whether each component matches.
        bool equal = true;
        for (int j = 0; equal && j < sz; j++)
          {
          equal &= this->PointSet[h][j] == temp[j];
          }
        // If matched, the polytope is present at index h.
        if (equal) 
          {
          break;
          }
        } 
      // No match: advance to next slot in table.
      h++;
      if (h == PolytopeGeometry::TABSIZE) h = 0;
      }
    }

  // Guaranteed to have found or created the polytope in slot h.
  return this->Index[h];
} // AddPtope

// ---------------------------------------------------------------------

// Return a domain coordinate for an active vertex via a
// caller-allocated buffer.
void PolytopeGeometry::GetActiveDomCoord(IdType id, double *c)
{
  // pass buffer to PT object to fill.
  this->ActiveDomCoord->GetPoint(-1-id, c);

}

// ------------------------------------------------------------------

// Copy a range coordinate into caller-allocated buffer.
// TODO: consider removing need to copy by passing reference to 
//       coordinate data in VTK array.
void PolytopeGeometry::GetActiveRngCoord(IdType id, double *c)
{
  id = -1 - id;
  for (int i = 0; i < this->N; i++)
    {
    c[i] = this->ActiveScalars[id][i];
    }
}


// ------------------------------------------------------------------

// Return a specified component of an active vertex range coordinate.
double PolytopeGeometry::GetActiveRngCoordComponent(IdType id,  int c)
{
  id = -1 -id; 
  return this->ActiveScalars[id][c];
}

// ------------------------------------------------------------------

double  PolytopeGeometry::GetStoredRangeComponent(IdType id, int c){
  return StoredScalars[id][c];
}

// ------------------------------------------------------------------

// Split an active polytope against a sequence of increasing 
// thresholds, returning new polytope fragments and plane-polytope 
// intersections ("ends") via caller-allocated buffers.  Note
// that number of returned components will be one more than number
// of planes, due to a possible residue beyond the last plane.
// "nrts" refers to the number of planes, and does not include this
// extra entry.

void PolytopeGeometry::Split(
  IdType ptid,           // polytope id to split
  double *thresholds,    // ascending sequence of thresholds.
  int nrts,              // number of thresholds in buffer
  int field,             // component of domain being sliced
  IdType *frags,         // return buffer for fragments
  IdType *ends)          // return buffer for intersections
{
  struct _Ptope *entry;
  int t_comp, prev_comp;
  IdType a;
  double ai;

  bool placed = false;    // has the ptope being fully processed?
 

  entry = this->Ptope + ptid;
  
  // Easy case: reuse cached result.
  if (entry->fragments >= 0 && entry->ends >= 0 && entry->splitDim == field)  
    {
    // This ptope has already been fragmented on this axis, 
    // so we re-use the fragment and end information cached with the ptope entry.
    for (int i = 0; i < nrts+1; i++)
      {
      frags[i] = this->RecCache[entry->fragments + i];
      ends[i] = this->RecCache[entry->ends + i];
      }
    return;
    }


  // Allocate space for fragments and ends in the global cache, enlarging
  // the cache if necessary. Combine frag/end space to reduce reallocation.
  // NB: allocate space for one fragment per cutting plane PLUS residue.

  this->EnsureCacheCapacity(this->CacheSize + 2*(nrts + 1));
  entry->fragments = this->CacheSize;
  entry->ends = entry->fragments + nrts + 1;
  entry->splitDim = field;
  this->CacheSize += 2*(nrts + 1);

  // Handle polytope based on dimension.
  switch (entry->dim) {

    case 0: 

      // Rare: polytope is a single point.
      // Point either:
      // - lies on a plane
      // - lies between two planes
      // - comes after the last plane.

      a = this->Facet[entry->start];
      ai = this->GetActiveRngCoordComponent(a, field);

      prev_comp = -1; // comparison result for previous plane. 
      for (int c = 0; c < nrts; c++)
        {
        // Compare coordinate component with current plane threshold.
        t_comp = Compare(ai, thresholds[c]);
        switch (t_comp) {
          case -1: 
            // plane is > component.
            this->RecCache[entry->fragments + c] = this->RecCache[entry->ends + c] = ends[c] = frags[c] = 0; 
            break;
          case  0: 
            // component matches plane, point lies on plane.
            this->RecCache[entry->ends + c] = ends[c] = a; 
            placed = true;
            this->RecCache[entry->fragments + c] = frags[c] = 0; 
            break;
          case  1:
            // component is < plane 
            if (prev_comp < 0)
              {
              // previous plane was < component, so point
              // lies between both planes.
              this->RecCache[entry->ends + c] = ends[c] = 0; 
              this->RecCache[entry->fragments + c] = frags[c] = a;
              placed = true;
              }
            else
              {
              this->RecCache[entry->ends + c] = this->RecCache[entry->fragments + c] = ends[c] = frags[c] = 0;
              }
            break;
          }
        prev_comp = t_comp;
        }
      // handle case when point lies beyond last cutting plane.
      this->RecCache[entry->ends + nrts] = ends[nrts] 
          = this->RecCache[entry->fragments + nrts] = frags[nrts]
          = placed ? 0 : a;
      // END OF ZERO-D (point) CASE
      break;
      
    case 1: // EDGE 
      this->SplitEdge(entry, thresholds, nrts, field, frags, ends);
      break;

    default: // arbitrary dimension polytope.
      this->SplitPtope(entry, thresholds, nrts, field, frags, ends);

    } // switch dimension

} // Split.




// ------------------------------------------------------------------

// Split an edge against a set of cutting planes to produce a 
// collection of segments.  Degenerate cases, e.g. where an edge
// lies within a plane or between planes, are handled.

void PolytopeGeometry::SplitEdge(
  struct _Ptope *entry,   // The polytope to be split (will be an edge)
  double *thresholds,     // Thresholds on which to cut the edge
  int nrts,               // Number of thresholds provided
  int field,              // Coordinate component being cut.
  IdType *frags,          // Buffer for fragments (includes a residue slot)
  IdType *ends)           // Buffer for intersections
{
  double t;
  int t_comp, prev_comp;
  IdType temp[2];      // temporary edge between two vertices.
  IdType a, b;

  // Buffers for domain/range coordinates of edge endpoints.
  double domCoordA[MaxDomainDim];
  double domCoordB[MaxDomainDim];

  double rngCoordA[MaxRangeDim];
  double rngCoordB[MaxRangeDim];

  // Buffers for interpolated coordinates.
  double domCoord[MaxDomainDim];
  double rngCoord[MaxRangeDim];
  
  double ai, bi;
  IdType coord_prev;
  bool working;

  // Get indices of edge endpoints, range coordinates of
  // endpoints, and values for component being cut against.
  a = this->Facet[entry->start];
  b = this->Facet[entry->start+1];
  this->GetActiveRngCoord(a, rngCoordA);
  this->GetActiveRngCoord(b, rngCoordB);
  ai = rngCoordA[field];
  bi = rngCoordB[field];

  // Ensure edge is oriented so that ai < bi
  if (bi < ai)
    {
    std::swap<double>(ai, bi);
    std::swap<IdType>(a, b);
    for (int i = 0; i < this->N; i++)
      {
      std::swap<double>(rngCoordA[i], rngCoordB[i]);
      }
    }

  // Separate out the case of edge lying in the cutting plane.
  
  if (Compare(ai, bi))
    { 
    // edge not coplanar, so seek intersections.  

    coord_prev = a;   // Last endpoint of edge to be processed
    working = true;   // Is output of edge ongoing?
    
    for (int c = 0; c < nrts; c++)
      {
      // Find parametric coordinate of edge/plane intersection,
      // and classify as below / on / above lowest endpoint.
      t = (thresholds[c] - ai) / (bi - ai);
      t_comp = Compare(t, 0.0);

      switch (t_comp) {
        case -1: 
          // edge starts at component value above that of the plane.
          this->RecCache[entry->fragments + c] = this->RecCache[entry->ends + c] = frags[c] = ends[c] = 0;
          break;

        case  0: 
          // left-end of edge lies on the plane
          this->RecCache[entry->fragments + c] = this->RecCache[entry->ends + c] = frags[c] = ends[c] = 0; 
         break;                

        case  1: 
          // edge starts below plane.
          // consider intersection with right hand of edge.
          t_comp = Compare(t, 1.0);             
          switch (t_comp) {
            case -1: 
              // Edge ends at component value above plane, so is cut.
              // Interpolate intersection points in both 
              // domain and range spaces.
              this->GetActiveDomCoord(a, domCoordA);
              this->GetActiveDomCoord(b, domCoordB);
           

              for (int j = 0; j < this->M_topo+1; j++)
                {
                domCoord[j] = (1-t)*domCoordA[j] + t*domCoordB[j];
                }
              for (int j = 0; j < this->N; j++)
                {
                rngCoord[j] = (1-t)*rngCoordA[j] + t*rngCoordB[j];
                }           
              // Edge segment starts at the previous edge endpoint computed.         
              temp[0]  = coord_prev;
              // Edge segment ends at intersection point.  This intersection
              // point is also returned through the end/intersection buffer.
              this->RecCache[entry->ends + c] = ends[c]  = temp[1] = this->AddCoord(domCoord, rngCoord);
              // New polytope is defined for the edge fragment.
              this->RecCache[entry->fragments + c] = frags[c] = this->AddPtope(temp, 1, 2);
              // This intersection point now becomes the last edge endpoint seen.
              coord_prev = temp[1];
              break;

            case  0: 
              // edge ends at cutting plane
              // connect previous edge endpoint to right end of edge.
              temp[0]  = coord_prev;
              temp[1]  = b;
              this->RecCache[entry->fragments + c] = frags[c] = this->AddPtope(temp, 1, 2);
              this->RecCache[entry->ends + c] = ends[c]  = b;
              coord_prev = temp[1];
              // as edge endpoint has been seen, there will be no residue.
              working = false;
              break;
              
            case  1: 
              // edge ends at a component value below that of the 
              // current plane. If the end of the edge has not been seen,
              // there is still a segment to output.

              if (working) 
                {
                // connect previous endpoint to endpoint of edge.
                temp[0]  = coord_prev;
                temp[1]  = b;
                this->RecCache[entry->fragments + c] = frags[c] = this->AddPtope(temp, 1, 2);
                working  = false;                    
                }
              else
                {
                this->RecCache[entry->fragments + c] = frags[c] = 0;
                }
              this->RecCache[entry->ends + c] = ends[c]  = 0;
              break;
            } // switch on right-end.

        } // switch on left-end.
        
      } // for each cutting plane

    // If we haven't reached the end of the line, there is a
    // residue beyond the last cutting plane.  
    if (working)
      {
      temp[0] = coord_prev;
      temp[1] = b;
      this->RecCache[entry->fragments + nrts] = frags[nrts] = this->AddPtope(temp, 1, 2);
      // NB since line end didn't hit a plane, end is empty.
      this->RecCache[entry->ends + nrts] = ends[nrts] = 0;
      }
    else
      {
      this->RecCache[entry->fragments + nrts] = this->RecCache[entry->ends + nrts] = frags[nrts] = ends[nrts] = 0;  
      }
      
    } // edge not parallel to planes
  else 
    { 
    // Edge lies in a plane parallel to the cutting planes for this component.
    // Scan through planes until the position of the edge is determined.
    prev_comp = -1;
    working = false;
    for (int c = 0; c < nrts; c++)
      {
      t_comp = Compare(thresholds[c], ai);
      switch (t_comp) {
        case -1: 
          this->RecCache[entry->ends + c] = this->RecCache[entry->fragments + c] = ends[c] = frags[c] = 0; 
          break;
        case  0: 
          this->RecCache[entry->ends + c] = ends[c] = (entry - this->Ptope);  // id of original ptope
          this->RecCache[entry->fragments + c] = frags[c] = 0; 
          working = true;
          break;
        case  1: 
          if (prev_comp < 0)
            {
            this->RecCache[entry->ends + c] = ends[c] = 0; 
            this->RecCache[entry->fragments + c] = frags[c] = (entry - this->Ptope);
            working = true;
            }
          else
            {
            this->RecCache[entry->ends + c] = this->RecCache[entry->fragments + c] = ends[c] = frags[c] = 0;
            }
          break;
        }
      prev_comp = t_comp;
      }
    // by definition, residue edges parallel to cutting plane does not 
    // intersect any cutting plane.
    this->RecCache[entry->ends + nrts] = ends[nrts] = 0;
    if (working)
      {
      this->RecCache[entry->fragments + nrts] = frags[nrts] = 0;
      }
    else
      {
      this->RecCache[entry->fragments + nrts] = frags[nrts] = (entry - this->Ptope); 
      }
    } // else edge parallel to cutting plane

} // SplitEdge.

// ------------------------------------------------------------------

// Split a polytope of dimension d > 1.
// Recursively split the polytope's facets, then
// at each threshold,
// - attempt to build a new "end" facet from the fragments
//   of the facet intersections
// - attempt to build a new fragment from:
//   1.  the fragments of the facets
//   2.  the end facet, if it exists
//   3.  the PREVIOUS end facet, if it exists.

void PolytopeGeometry::SplitPtope(
  struct _Ptope *entry,   // the ptope to be split
  double *thresholds,     // cutting planes, ascending
  int nrts,               // number of cutting planes
  int field,              // component being cut
  IdType *frags,       // buffer to hold fragments
  IdType *ends)        // buffer to hold intersections.
{

  // Allocate temp space for storing output fragments and ends formed by 
  // cutting each part of this ptope against the thresholds.
  // There are nrts + 1 possible sets of output, as some parts of the
  // ptope may lie beyond the final threshold.
  // There are entrysize + 2 rows, as in addition to fragments lying 
  // between cutting planes, there are up to two end ptopes to consider.

  IdType nFragNr= (nrts + 1) * (entry->size + 2);
/*  this->EnsureCacheCapacity(this->CacheSize + 2*nFragNr);
  IdType *rec_frags = this->RecCache + this->CacheSize;
  IdType *rec_ends  = rec_frags + nFragNr; 
  this->CacheSize += 2*nFragNr;*/
  
  IdType * rec_frags = new IdType [nFragNr]; 
  IdType * rec_ends = new IdType [nFragNr]; 



  // ids for the polytopes formed by fragments lying in cutting planes.
  IdType end_tope, end_prev;
  int nrParts;
  
  // Recursively split each facet on the same field/thresholds.
  // Note offsets to start of buffer for each threshold include
  // the slots for the residues.
  for (int f = 0; f < entry->size; f++)
    {
    this->Split(
        this->Facet[entry->start + f],
        thresholds,      
        nrts,
        field,
        rec_frags + f*(nrts + 1),
        rec_ends  + f*(nrts + 1));
    } // for each face
  

  // Try to compute a new polytope formed by each cutting 
  // plane and by the residue (hence <= condition for loop 
  // termination).
  // We convert end fragments into a new polytope, adding
  // these to each pair of polytopes, using the next
  // two rows of the rec_ table ... so assumption is that
  // input polytopes are at least to facets smaller than
  // MaxNrFacetsPerPtope.

  end_prev = 0;
  for (int i = 0; i <= nrts; i++)
    {
    // Track number of *candidate* parts for each new polytope.
    // Each facet will contribute a candidate.
    nrParts = entry->size;    
    // Attempt to build an intersection polytope from end fragments
    // lying in the plane (for the residue this will always be nil).
    end_tope = this->RecCache[entry->ends + i] = ends[i]   
       = this->Select(rec_ends, nrts+1, nrParts, i, entry->dim-1, entry->dim);
    // If there is an intersection polytope, add it to the
    // candidate components for making the next polytope fragment.
    if (end_tope)
      {
      rec_frags[(nrts + 1) * nrParts + i] = end_tope;
      nrParts++;
      }
    // If the *previous* cut resulted in an end-polytope, add it
    // to the candidate fragments for this polytope.
    if (end_prev)
      {
      rec_frags[(nrts + 1) * nrParts + i] = end_prev;  
      nrParts++;
      }
    // The end polytope from this plane will be the previous for the next plane.
    end_prev = end_tope;
    // Attempt to build a polytope fragment from the facet fragments and any 
    // end polytopes.
    this->RecCache[entry->fragments + i] = frags[i] 
        = this->Select(rec_frags, nrts+1, nrParts, i, entry->dim, entry->dim);
    }
    
    delete rec_frags;
    delete rec_ends;

    
} // SplitPtope


// ------------------------------------------------------------------

// Identify a polytope of dimension between d and endDim, where
// endDim <= d, from a collection of polytope fragments held in 
// a 2D array pointed to by "pts", with "nrts" entries per row.
// Polytope fragments are to be taken from an array column 
// given by "slice".  Rationale for this data layout is to be
// be found in the implementation of "Split" given below.
// NOTE: the "nrts" parameter passed  to Select INCLUDES the
// extra slot allocated for the residue (see Split). 

IdType PolytopeGeometry::Select(
  IdType *pts,
  int nrts, // NB here = number of planes + 1 for the residue
  int nrps,
  int slice,
  int d, int endDim)
{  
  struct _Ptope *entry;
  IdType pt;
  int next, pos;
  IdType temp[MaxNrSlicesPerPtope];
  int end = endDim > 0 ? endDim : 0;

  // Attempt to form the largest possible face.

  // For efficiency, separate cases for searching for points or
  // searching for higher-dim polytopes.
  // d = dimensionality of components sought.
  
  // EITHER we find a single ptope of the required dimension,
  // OR we find sufficient facets to build one.

  // 1. Attempt to find one ptope:
  next = 0;
  for (int i = 0; i < nrps; i++)
    {
    if ((pt = pts[i*nrts + slice]) >= 0 && this->Ptope[pt].dim == d)
      {
      temp[next++] = pt;
      }
    } // for each candidate component
  
  if (next == 1) { return temp[0]; }  // exactly one, existing, ptope found.
  
  // 2. Attempt to build from facets.
  // Separate case of 1D (lines from points) & 2D; as vertices are stored
  // as negative values their dimension is implicit, saving on inner-loop 
  // tests.


  next = 0;
  if (d > 1)
    { 
    // not a line, so must test candidate dimension explicitly.
    for (int i = 0; i < nrps; i++)
      {
      if ((pt = pts[i*nrts + slice]) > 0 && this->Ptope[pt].dim == d-1)
        {
        for (pos = 0; pos < next && temp[pos] != pt; pos++) {;} 
        if (pos == next)
          { 
          temp[next++] = pt;
          }
        }
      } // for each candidate point
    if (next > d) 
      { 
      // create a new ptope of one dimension up from its component faces.
      return this->AddPtope(temp, d, next); 
      }
    } // if not a line.
  else
    { 
    // working with a 1-D line, look for -ve components.
    for (int i = 0; i < nrps; i++)
      {
      if ((pt = pts[i*nrts + slice]) < 0)
        {
        for (pos = 0; pos < next && temp[pos] != pt; pos++) {;} 
        if (pos == next)
          { 
          temp[next++] = pt;
          }
        }
      } // for each candidate point
    if (next == 2)
      {
      return this->AddPtope(temp, 1, 2); 
      }
    } // line case.

  // The worst case: we have no suitable polytope;
  // return the degenerate Nil polytope.
  return 0;
}



// ------------------------------------------------------------------

// Return the pre-computed lower corner of the bounding box for
// an active polytope.
double *PolytopeGeometry::GetPtopeCoordMin(IdType id)
{
  return this->DblCache + this->Ptope[id].minCoord;
}

// ------------------------------------------------------------------

// ------------------------------------------------------------------

// Return pointer to buffer holding the center of a given polytope;
// center is cached in the polytope, and is computed on demand. 
// Optimize special case of low-dimensional polytopes, but in general
// invoke function recursively on facets, averaging over facet centers.
double *PolytopeGeometry::GetPtopeCenter(IdType pid)
{
  struct _Ptope *entry = this->Ptope + pid;
  double *fcenter, *base;
  double point[MaxDomainDim];
  IdType a;
  IdType *facet;

  if (entry->center >= 0)
    {
    // pre-computed point, return pointer to cache.
    return this->DblCache + entry->center;
    }   
    
  this->EnsureDblCacheCapacity(this->DblCacheSize + this->M_geom);
  // Changed to multi-dimensional point
  entry->center = this->DblCacheSize;
  this->DblCacheSize += this->M_geom;
    
  // No cached center, so compute cache.
  // Initialize cache entry to zeroes.
  for (int i = 0; i < this->M_geom; i++)
    {    
    this->DblCache[entry->center + i] = 0.0;
    }

    
  // Branch on polytope dimension to handle common cases cheaply.    
  facet = this->Facet + this->Ptope[pid].start;
  switch (entry->dim) {

    case 0:
      // Point: center is just the point itself.
      a = -1 - facet[0];
      this->ActiveDomCoord->GetPoint(a, this->DblCache + entry->center);
      break;
      
    case 1:
      // Edge: center is midway between endpoints.
      a = -1 - facet[0];
      this->ActiveDomCoord->GetPoint(a, point);
      a = -1 - facet[1];
      this->ActiveDomCoord->GetPoint(a, this->DblCache + entry->center);
 
      for (int j = 0; j < this->M_geom; j++)
        {
        this->DblCache[entry->center + j] = 0.5 * (this->DblCache[entry->center + j] + point[j]);
        }

      break;
      
    default:
      // Recursively compute facet centers and 
      // accumulate sum.
         
      for (int i = 0; i < entry->size; i++)
        {
        fcenter = this->GetPtopeCenter(facet[i]);
        for (int j = 0; j < this->M_geom; j++)
          {
          this->DblCache[entry->center + j] += fcenter[j];
          }
        }
      
      // Take average of facet centers
      // as polytope center.
      for (int j = 0; j < this->M_geom; j++)
        {
        this->DblCache[entry->center + j] /= entry->size;
        }
    }
  return this->DblCache + entry->center;
}


// ------------------------------------------------------------------

void PolytopeGeometry::EnsureCacheCapacity(IdType size)
{
  if (this->MaxCache >= size)
    return;

  IdType *newCache;
  IdType newSize = this->MaxCache ? this->MaxCache*2 : 8192;

  newCache = new IdType [newSize];
  if (!newCache){
    MsgAndExit("Cannot allocate polytope cache.");
    return;
  }
  std::memcpy(newCache, RecCache, (this->CacheSize)*(sizeof (IdType)));
  
  delete [] this->RecCache;
  this->RecCache = newCache;
  this->MaxCache = newSize;
}


// ------------------------------------------------------------------

void PolytopeGeometry::EnsureDblCacheCapacity(IdType size)
{
  if (this->MaxDblCache >= size)
    return;

  double *newCache;
  IdType newSize = this->MaxDblCache ? this->MaxDblCache*2 : 8192;

  newCache = new double [newSize];
  if (!newCache)  {
    MsgAndExit("Cannot allocate double cache.");
    return;
  }
  
  std::memcpy(newCache, DblCache, (this->DblCacheSize)*(sizeof (double))); 

  delete [] this->DblCache;
  this->DblCache = newCache;
  this->MaxDblCache = newSize;
}
