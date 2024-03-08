#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Part 1.
      std::vector<Vector2D> newPoints;
      //newPoints.push_back(points.front());
      for (int i = 0; i < points.size() - 1; i++) {
          Vector2D newPoint = t * points[i] + (1 - t) * points[i + 1];
          newPoints.push_back(newPoint);
      }
      //newPoints.push_back(points.back());
      return newPoints;
      //return std::vector<Vector2D>();
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
      std::vector<Vector3D> newPoints;
      //newPoints.push_back(points.front());
      for (int i = 0; i < points.size() - 1; i++) {
          Vector3D newPoint = t * points[i] + (1 - t) * points[i + 1];
          newPoints.push_back(newPoint);
      }
      //newPoints.push_back(points.back());
      return newPoints;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
      std::vector<Vector3D> intermediate = evaluateStep(points, t);
      //double btw_point = (intermediate.size() - 1) * t;
      //double t1 = btw_point - floor(btw_point);
      //return intermediate[floor(btw_point)];
      /*
      if (btw_point == intermediate.size() - 1) {
          return intermediate[floor(btw_point)];
      }
      else {
          return t1 * intermediate[floor(btw_point)] + (1 - t1) * intermediate[floor(btw_point) + 1];
      }
      */
      //return t1 * intermediate[floor(btw_point)] + (1 - t1) * intermediate[floor(btw_point) + 1];
      while (intermediate.size() > 1) {
          intermediate = evaluateStep(intermediate, t);
      }
      return intermediate[0]; 
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Part 2.
      std::vector<Vector3D> v_curve;
      for (int i = 0; i < controlPoints.size(); i++) {
          Vector3D u_point = evaluate1D(controlPoints[i], u);
          v_curve.push_back(u_point);
      }
      Vector3D uv_point = evaluate1D(v_curve, v);

      return uv_point;
  }

  Vector3D face_normal(HalfedgeCIter h) {
      //only supports triangles
      VertexCIter v0 = h->vertex();
      VertexCIter v1 = h->next()->vertex();
      VertexCIter v2 = h->next()->next()->vertex();
      Vector3D v01 = v1->position - v0->position;
      Vector3D v12 = v2->position - v1->position;
      return 0.5 * cross(v01, v12);
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
      Vector3D sum_vec;
      HalfedgeCIter h = _halfedge;      // get the outgoing half-edge of the vertex
      do {
          sum_vec += face_normal(h);
          HalfedgeCIter h_twin = h->twin(); // get the opposite half-edge
          //VertexCIter v = h_twin->vertex(); // vertex is the 'source' of the half-edge, so
                                              // h->vertex() is v, whereas h_twin->vertex()
                                              // is the neighboring vertex
          //cout << v->position << endl;      // print the vertex position
          h = h_twin->next();               // move to the next outgoing half-edge of the vertex

      } while (h != this->halfedge());          // keep going until we are back where we were
      return sum_vec / sum_vec.norm() + position;
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
      std::cout << "test!!" << std::endl;
      if (e0->isBoundary()) {
          // Can't flip a boundary edge
          std::cout << "can't flip" << std::endl;
          return e0;
      }
      HalfedgeIter h0 = e0->halfedge();
      HalfedgeIter h1 = h0->next();
      HalfedgeIter h2 = h1->next();
      HalfedgeIter h3 = h0->twin();
      HalfedgeIter h4 = h3->next();
      HalfedgeIter h5 = h4->next();
      HalfedgeIter h6 = h1->twin();
      HalfedgeIter h7 = h2->twin();
      HalfedgeIter h8 = h4->twin();
      HalfedgeIter h9 = h5->twin();

      VertexIter v0 = h0->vertex();
      VertexIter v1 = h3->vertex();
      VertexIter v2 = h2->vertex();
      VertexIter v3 = h5->vertex();

      EdgeIter e1 = h1->edge();
      EdgeIter e2 = h2->edge();
      EdgeIter e3 = h4->edge();
      EdgeIter e4 = h5->edge();

      FaceIter f0 = h0->face();
      FaceIter f1 = h3->face();


      check_for(h1);

      ///*
      

      h0->next() = h1;
      h0->twin() = h3;
      h0->vertex() = v3;
      h0->edge() = e0;
      h0->face() = f0;
      h1->next() = h2;
      h1->twin() = h7;
      h1->vertex() = v2;
      h1->edge() = e2;
      h1->face() = f0;
      h2->next() = h0;
      h2->twin() = h8;
      h2->vertex() = v0;
      h2->edge() = e3;
      h2->face() = f0;
      h3->next() = h4;
      h3->twin() = h0;
      h3->vertex() = v2;
      h3->edge() = e0;
      h3->face() = f1;
      h4->next() = h5;
      h4->twin() = h9;
      h4->vertex() = v3;
      h4->edge() = e4;
      h4->face() = f1;
      h5->next() = h3;
      h5->twin() = h6;
      h5->vertex() = v1;
      h5->edge() = e1;
      h5->face() = f1;
      h6->next() = h6->next();
      h6->twin() = h5;
      h6->vertex() = v2;
      h6->edge() = e1;
      h6->face() = h6->face();
      h7->next() = h7->next();
      h7->twin() = h1;
      h7->vertex() = v0;
      h7->edge() = e2;
      h7->face() = h7->face();
      h8->next() = h8->next();
      h8->twin() = h2;
      h8->vertex() = v3;
      h8->edge() = e3;
      h8->face() = h8->face();
      h9->next() = h9->next(); // didn’t change, but set it anyway!
      h9->twin() = h4;
      h9->vertex() = v1;
      h9->edge() = e4;
      h9->face() = h9->face(); // didn’t change, but set it anyway!
      
      
      v0->halfedge() = h2;
      v1->halfedge() = h5;
      v2->halfedge() = h3;
      v3->halfedge() = h0;
      
      e0->halfedge() = h0;
      e1->halfedge() = h5;
      e2->halfedge() = h1;
      e3->halfedge() = h2;
      e4->halfedge() = h4;

      f0->halfedge() = h0;
      f1->halfedge() = h3;

      
      //*/
      return e0;
    
      
  }




  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.

      if (e0->isBoundary()) {
          // Can't flip a boundary edge
          std::cout << "can't flip" << std::endl;
          return VertexIter();
      }

      

      HalfedgeIter h0 = e0->halfedge();
      HalfedgeIter h1 = h0->next();
      HalfedgeIter h2 = h1->next();
      HalfedgeIter h3 = h0->twin();
      HalfedgeIter h4 = h3->next();
      HalfedgeIter h5 = h4->next();
      HalfedgeIter h6 = h1->twin();
      HalfedgeIter h7 = h2->twin();
      HalfedgeIter h8 = h4->twin();
      HalfedgeIter h9 = h5->twin();
      HalfedgeIter h10 = newHalfedge();
      
      HalfedgeIter h11 = newHalfedge();
      
      HalfedgeIter h12 = newHalfedge();
      
      HalfedgeIter h13 = newHalfedge();
      
      HalfedgeIter h14 = newHalfedge();
      
      HalfedgeIter h15 = newHalfedge();
      

      // check_for(h1);


      VertexIter v0 = h0->vertex();
      VertexIter v1 = h3->vertex();
      VertexIter v2 = h2->vertex();
      VertexIter v3 = h5->vertex();
      VertexIter v4 = newVertex();
      v4->halfedge() = h0;

      e0->isNew = false();
      EdgeIter e1 = h1->edge();
      EdgeIter e2 = h2->edge();
      EdgeIter e3 = h4->edge();
      EdgeIter e4 = h5->edge();
      EdgeIter e5 = newEdge();
      e5->halfedge() = h0;
      e5->isNew = false;
      EdgeIter e6 = newEdge();
      e6->halfedge() = h0;
      e5->isNew = true;
      EdgeIter e7 = newEdge();
      e7->halfedge() = h0;
      e5->isNew = true;
      EdgeIter e8 = newEdge();
      e8->halfedge() = h0;

      FaceIter f0 = h0->face();
      FaceIter f1 = h3->face();
      FaceIter f2 = newFace();
      f2->halfedge() = h0;
      FaceIter f3 = newFace();
      f3->halfedge() = h0;

     // bool b = checkFaceValidity(f0);

      
      
      



      ///*


      h0->next() = h14;
      h0->twin() = h3;
      h0->vertex() = v0;
      h0->edge() = e8;
      h0->face() = f0;
      h1->next() = h0;
      h1->twin() = h7;
      h1->vertex() = v2;
      h1->edge() = e2;
      h1->face() = f0;
      h2->next() = h13;
      h2->twin() = h8;
      h2->vertex() = v0;
      h2->edge() = e3;
      h2->face() = f1;
      h3->next() = h2;
      h3->twin() = h0;
      h3->vertex() = v4;
      h3->edge() = e8;
      h3->face() = f1;
      h4->next() = h11;
      h4->twin() = h9;
      h4->vertex() = v3;
      h4->edge() = e4;
      h4->face() = f2;
      h5->next() = h15;
      h5->twin() = h6;
      h5->vertex() = v1;
      h5->edge() = e1;
      h5->face() = f3;
      h6->next() = h6->next();
      h6->twin() = h5;
      h6->vertex() = v2;
      h6->edge() = e1;
      h6->face() = h6->face();
      h7->next() = h7->next();
      h7->twin() = h1;
      h7->vertex() = v0;
      h7->edge() = e2;
      h7->face() = h7->face();
      h8->next() = h8->next();
      h8->twin() = h2;
      h8->vertex() = v3;
      h8->edge() = e3;
      h8->face() = h8->face();
      h9->next() = h9->next(); // didn’t change, but set it anyway!
      h9->twin() = h4;
      h9->vertex() = v1;
      h9->edge() = e4;
      h9->face() = h9->face(); // didn’t change, but set it anyway!
      h10->next() = h5;
      h10->twin() = h11;
      h10->vertex() = v4;
      h10->edge() = e5;
      h10->face() = f3;
      h11->next() = h12;
      h11->twin() = h10;
      h11->vertex() = v1;
      h11->edge() = e5;
      h11->face() = f2;
      h12->next() = h4;
      h12->twin() = h13;
      h12->vertex() = v4;
      h12->edge() = e6;
      h12->face() = f2;
      h13->next() = h3;
      h13->twin() = h12;
      h13->vertex() = v3;
      h13->edge() = e6; 
      h13->face() = f1;
      h14->next() = h1;
      h14->twin() = h15;
      h14->vertex() = v4;
      h14->edge() = e7;
      h14->face() = f0;
      h15->next() = h10;
      h15->twin() = h14;
      h15->vertex() = v2;
      h15->edge() = e7;
      h15->face() = f3;
      v0->halfedge() = h2;
      v1->halfedge() = h5;
      v2->halfedge() = h1;
      v3->halfedge() = h4;

      v4->halfedge() = h3;

      e0->halfedge() = h0;
      e1->halfedge() = h5;
      e2->halfedge() = h1;
      e3->halfedge() = h2;
      e4->halfedge() = h4;
      e5->halfedge() = h10;
      e6->halfedge() = h12;
      e7->halfedge() = h14;
      e8->halfedge() = h0;

      f0->halfedge() = h0;
      f1->halfedge() = h13;
      f2->halfedge() = h11;
      f3->halfedge() = h15;

      //deleteEdge(e0);
      Vector3D pos = (v0->position + v1->position) / 2;
      v4->position = pos;

      //checkFaceValidity(f0);

      //*/
      return v4;

  }

  /*
  
  bool checkFaceValidity(FaceIter& f) {
      if (f->isBoundary()) {
          //std::cerr << "Face is a boundary, which does not have a closed halfedge loop.\n";
          return false;
      }

      HalfedgeIter start = f->halfedge(); // Starting halfedge of the face
      if (!start) {
          return false;
      }

      // Use a set to track visited halfedges to detect loops correctly
      std::set<Halfedge*> visited;
      HalfedgeIter h = start;
      do {
          // Check if we've visited this halfedge already (which would indicate a loop)
          if (visited.find(&*h) != visited.end()) {
              return false;
          }
          visited.insert(&*h);

          // Check if the halfedge correctly points back to the face
          if (h->face() != f) {
              return false;
          }

          // Check if the twin's twin is the halfedge itself
          if (h->twin()->twin() != h) {
              return false;
          }

          // Check if next halfedge is valid (not null)
          if (!h->next()) {
              return false;
          }

          // Move to the next halfedge in the face
          h = h->next();
      } while (h != start); // Loop until we return to the starting halfedge

      // If all checks passed, the face and its associated elements are valid
      return true;
  }


  */

  


  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.

      /*
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          if (e->isNew == false) {
              VertexIter v = mesh.splitEdge(e);
              e->isNew = true;
          }

      }
      */

      

      ///*
      
      std::cout << "hihi" << std::endl;
      for (FaceIter f = mesh.facesBegin(); f != mesh.facesEnd(); f++) {
      //FaceIter f = mesh.facesBegin();
          HalfedgeIter h0 = f->halfedge();
          HalfedgeIter h1 = h0->next();
          HalfedgeIter h2 = h1->next();

          VertexIter v0 = h0->vertex();
          VertexIter v1 = h1->vertex();
          VertexIter v2 = h2->vertex();

          v0->isNew = false;
          v1->isNew = false;
          v2->isNew = false;

          Vector3D pos = (v2->position + v1->position) / 2;
          v0->newPosition = pos;
          //h1->edge->newPosition = pos;
          EdgeIter e = h1->edge();
          if (e->isNew == false) {
              VertexIter v = mesh.splitEdge(e);
              e->isNew = true;
          }
      }

      
      //*/

      
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    
    // 4. Flip any new edge that connects an old and new vertex.

    // 5. Copy the new vertex positions into final Vertex::position.
      
  }
}
