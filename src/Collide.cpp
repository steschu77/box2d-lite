/*
 * Copyright (c) 2006-2007 Erin Catto http://www.gphysics.com
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies.
 * Erin Catto makes no representations about the suitability
 * of this software for any purpose.
 * It is provided "as is" without express or implied warranty.
 */

#include "Arbiter.h"
#include "Body.h"

// Box vertex and edge numbering:
//
//        ^ y
//        |
//        e1
//   v2 ------ v1
//    |        |
// e2 |        | e4  --> x
//    |        |
//   v3 ------ v4
//        e3

enum Axis
{
  FACE_A_X,
  FACE_A_Y,
  FACE_B_X,
  FACE_B_Y
};

enum EdgeNumbers
{
  NO_EDGE = 0,
  EDGE1,
  EDGE2,
  EDGE3,
  EDGE4
};

struct ClipVertex
{
  ClipVertex() { fp.value = 0; }
  x3d::vector2 v;
  FeaturePair fp;
};

template <typename T>
inline void Swap(T& a, T& b)
{
  T tmp = a;
  a = b;
  b = tmp;
}

void Flip(FeaturePair& fp)
{
  Swap(fp.e.inEdge1, fp.e.inEdge2);
  Swap(fp.e.outEdge1, fp.e.outEdge2);
}

struct ReferenceEdge
{
  ReferenceEdge(const Body* poly1, const Body* poly2, int flip);

  const Body* poly1;
  const Body* poly2;
  int flip;
  float separation;
  int index;
};

ReferenceEdge::ReferenceEdge(const Body* poly1, const Body* poly2, int flip)
: poly1(poly1)
, poly2(poly2)
, flip(flip)
, separation(-FLT_MAX)
, index(0)
{
}

// Find the max separation between poly1 and poly2 using edge normals from poly1.
static void b2FindMaxSeparation(ReferenceEdge* edge)
{
  const Body* poly1 = edge->poly1;
  const Body* poly2 = edge->poly2;
  const int count1 = poly1->count;
  const int count2 = poly2->count;
  const x3d::vector2* n1s = poly1->normals;
  const x3d::vector2* v1s = poly1->vertices;
  const x3d::vector2* v2s = poly2->vertices;

  const x3d::rot2 q = x3d::mulT(poly2->q, poly1->q);
  const x3d::vector2 p = x3d::mulT(poly2->q, poly1->p - poly2->p);

  for (int i = 0; i < count1; ++i) {

    // Get poly1 normal and vertex in poly2 coords.
    const x3d::vector2 n = x3d::mul(q, n1s[i]);
    const x3d::vector2 v1 = x3d::mul(q, v1s[i]) + p;

    // Find deepest point for normal i.
    float si = FLT_MAX;
    for (int j = 0; j < count2; ++j) {
      float sij = n * (v2s[j] - v1);
      if (sij < si) {
        si = sij;
      }
    }

    if (si > edge->separation) {
      edge->separation = si;
      edge->index = i;
    }
  }
}

static void b2FindIncidentEdge(ClipVertex c[2], ReferenceEdge* edge)
{
  const Body* poly1 = edge->poly1;
  const Body* poly2 = edge->poly2;

  const int count2 = poly2->count;
  const x3d::vector2* n2s = poly2->normals;
  const x3d::vector2* v2s = poly2->vertices;

  // Get the normal of the reference edge in poly2's frame.
  const x3d::vector2& n = poly1->normals[edge->index];
  const x3d::vector2 n1 = x3d::mulT(poly2->q, x3d::mul(poly1->q, n));

  // Find the incident edge on poly2.
  int index = 0;
  float minDot = FLT_MAX;
  for (int i = 0; i < count2; ++i) {
    float dot = n1 * n2s[i];
    if (dot < minDot) {
      minDot = dot;
      index = i;
    }
  }

  // Build the clip vertices for the incident edge.
  const int i1 = index;
  const int i2 = i1 + 1 < count2 ? i1 + 1 : 0;

  c[0].v = x3d::mul(poly2->q, v2s[i1]) + poly2->p;
  c[0].fp.e.inEdge2 = edge->index;
  c[0].fp.e.outEdge2 = i1;

  c[1].v = x3d::mul(poly2->q, v2s[i2]) + poly2->p;
  c[1].fp.e.inEdge2 = edge->index;
  c[1].fp.e.outEdge2 = i2;
}

int ClipSegmentToLine(ClipVertex vOut[2], ClipVertex vIn[2],
  const x3d::vector2& normal, float offset, char clipEdge)
{
  // Start with no output points
  int numOut = 0;

  // Calculate the distance of end points to the line
  float distance0 = normal * vIn[0].v - offset;
  float distance1 = normal * vIn[1].v - offset;

  // If the points are behind the plane
  if (distance0 <= 0.0f)
    vOut[numOut++] = vIn[0];
  if (distance1 <= 0.0f)
    vOut[numOut++] = vIn[1];

  // If the points are on different sides of the plane
  if (distance0 * distance1 < 0.0f) {
    // Find intersection point of edge and plane
    float interp = distance0 / (distance0 - distance1);
    vOut[numOut].v = vIn[0].v + interp * (vIn[1].v - vIn[0].v);
    if (distance0 > 0.0f) {
      vOut[numOut].fp = vIn[0].fp;
      vOut[numOut].fp.e.inEdge1 = clipEdge;
      vOut[numOut].fp.e.inEdge2 = NO_EDGE;
    } else {
      vOut[numOut].fp = vIn[1].fp;
      vOut[numOut].fp.e.outEdge1 = clipEdge;
      vOut[numOut].fp.e.outEdge2 = NO_EDGE;
    }
    ++numOut;
  }

  return numOut;
}

static void ComputeIncidentEdge(ClipVertex c[2], const x3d::vector2& h,
  const x3d::vector2& pos, const x3d::matrix2x2& Rot, const x3d::vector2& normal)
{
  // The normal is from the reference box. Convert it
  // to the incident boxe's frame and flip sign.
  x3d::matrix2x2 RotT = Rot.transpose();
  x3d::vector2 n = -(RotT * normal);
  x3d::vector2 nAbs = n.abs();

  if (nAbs.u0() > nAbs.u1()) {
    if (n.u0() > 0.0f) {
      c[0].v = x3d::vector2(h.u0(), -h.u1());
      c[0].fp.e.inEdge2 = EDGE3;
      c[0].fp.e.outEdge2 = EDGE4;

      c[1].v = x3d::vector2(h.u0(), h.u1());
      c[1].fp.e.inEdge2 = EDGE4;
      c[1].fp.e.outEdge2 = EDGE1;
    } else {
      c[0].v = x3d::vector2(-h.u0(), h.u1());
      c[0].fp.e.inEdge2 = EDGE1;
      c[0].fp.e.outEdge2 = EDGE2;

      c[1].v = x3d::vector2(-h.u0(), -h.u1());
      c[1].fp.e.inEdge2 = EDGE2;
      c[1].fp.e.outEdge2 = EDGE3;
    }
  } else {
    if (n.u1() > 0.0f) {
      c[0].v = x3d::vector2(h.u0(), h.u1());
      c[0].fp.e.inEdge2 = EDGE4;
      c[0].fp.e.outEdge2 = EDGE1;

      c[1].v = x3d::vector2(-h.u0(), h.u1());
      c[1].fp.e.inEdge2 = EDGE1;
      c[1].fp.e.outEdge2 = EDGE2;
    } else {
      c[0].v = x3d::vector2(-h.u0(), -h.u1());
      c[0].fp.e.inEdge2 = EDGE2;
      c[0].fp.e.outEdge2 = EDGE3;

      c[1].v = x3d::vector2(h.u0(), -h.u1());
      c[1].fp.e.inEdge2 = EDGE3;
      c[1].fp.e.outEdge2 = EDGE4;
    }
  }

  c[0].v = pos + Rot * c[0].v;
  c[1].v = pos + Rot * c[1].v;
}

// ----------------------------------------------------------------------------
template <typename T, size_t size>
size_t idxOfMin(const T (&x)[size])
{
  size_t idx = 0;

  for (size_t i = 1; i < size; i++) {
    if (x[i] < x[idx]) {
      idx = i;
    }
  }

  return idx;
}

// The normal points from A to B
int Collide(Contact* contacts, Body* bodyA, Body* bodyB)
{
  ReferenceEdge edgeA(bodyA, bodyB, 0);
  b2FindMaxSeparation(&edgeA);
  if (edgeA.separation > 0.0f) {
    return 0;
  }

  ReferenceEdge edgeB(bodyB, bodyA, 1);
  b2FindMaxSeparation(&edgeB);
  if (edgeB.separation > 0.0f) {
    return 0;
  }

  ReferenceEdge* edge = edgeB.separation > edgeA.separation ? &edgeB : &edgeA;

  ClipVertex incedge[2];
  b2FindIncidentEdge(incedge, edge);

	const int count1 = edge->poly1->count;
  const x3d::vector2* v1s = edge->poly1->vertices;

  const int iv1 = edge->index;
  const int iv2 = iv1 + 1 < count1 ? iv1 + 1 : 0;

  const x3d::vector2 v11 = x3d::mul(edge->poly1->q, v1s[iv1]) + edge->poly1->p;
  const x3d::vector2 v12 = x3d::mul(edge->poly1->q, v1s[iv2]) + edge->poly1->p;

	const x3d::vector2 tangent = (v12 - v11).norm();
  const x3d::vector2 xnormal = tangent.perpendicular();
  const float frontOffset = xnormal * v11;
  const float sideOffset1 = -tangent * v11 + 0.02f;
  const float sideOffset2 = tangent * v12 + 0.02f;

	ClipVertex incedge1[2];
  int np1 = ClipSegmentToLine(incedge1, incedge, -tangent, sideOffset1, iv1);

	ClipVertex incedge2[2];
  int np2 = ClipSegmentToLine(incedge2, incedge1, tangent, sideOffset2, iv2);

  // Setup
  x3d::vector2 hA = 0.5f * bodyA->width;
  x3d::vector2 hB = 0.5f * bodyB->width;

  x3d::vector2 posA = bodyA->position;
  x3d::vector2 posB = bodyB->position;

  x3d::matrix2x2 RotA(bodyA->rotation);
  x3d::matrix2x2 RotB(bodyB->rotation);

  x3d::matrix2x2 RotAT = RotA.transpose();
  x3d::matrix2x2 RotBT = RotB.transpose();

  x3d::vector2 dp = posB - posA;
  x3d::vector2 dA = RotAT * dp;
  x3d::vector2 dB = RotBT * dp;

  x3d::matrix2x2 C = RotAT * RotB;
  x3d::matrix2x2 absC = C.abs();
  x3d::matrix2x2 absCT = absC.transpose();

  // Box A faces
  x3d::vector2 faceA = dA.abs() - hA - absC * hB;
  x3d::vector2 faceB = dB.abs() - absCT * hA - hB;
  if (faceA.u0() > 0.0f || faceA.u1() > 0.0f || faceB.u0() > 0.0f
    || faceB.u1() > 0.0f) {
    return 0;
  }

  // Find best axis
  float separation[4];
  separation[0] = -faceA.u0();
  separation[1] = -faceA.u1();
  separation[2] = -faceB.u0();
  separation[3] = -faceB.u1();

  size_t idx = idxOfMin(separation);

  Axis axis = (Axis)idx;
  x3d::vector2 normal;

  // Setup clipping plane data based on the separating axis
  x3d::vector2 frontNormal, sideNormal;
  ClipVertex incidentEdge[2];
  float front, negSide, posSide;
  char negEdge, posEdge;

  // Compute the clipping lines and the line segment to be clipped.
  switch (axis) {
  case FACE_A_X: {
    normal = dA.u0() > 0.0f ? RotA.col(0) : -RotA.col(0);
    frontNormal = normal;
    front = posA * frontNormal + hA.u0();
    sideNormal = RotA.col(1);
    float side = posA * sideNormal;
    negSide = -side + hA.u1();
    posSide = side + hA.u1();
    negEdge = EDGE3;
    posEdge = EDGE1;
    ComputeIncidentEdge(incidentEdge, hB, posB, RotB, frontNormal);
  } break;

  case FACE_A_Y: {
    normal = dA.u1() > 0.0f ? RotA.col(1) : -RotA.col(1);
    frontNormal = normal;
    front = posA * frontNormal + hA.u1();
    sideNormal = RotA.col(0);
    float side = posA * sideNormal;
    negSide = -side + hA.u0();
    posSide = side + hA.u0();
    negEdge = EDGE2;
    posEdge = EDGE4;
    ComputeIncidentEdge(incidentEdge, hB, posB, RotB, frontNormal);
  } break;

  case FACE_B_X: {
    normal = dB.u0() > 0.0f ? RotB.col(0) : -RotB.col(0);
    frontNormal = -normal;
    front = posB * frontNormal + hB.u0();
    sideNormal = RotB.col(1);
    float side = posB * sideNormal;
    negSide = -side + hB.u1();
    posSide = side + hB.u1();
    negEdge = EDGE3;
    posEdge = EDGE1;
    ComputeIncidentEdge(incidentEdge, hA, posA, RotA, frontNormal);
  } break;

  case FACE_B_Y: {
    normal = dB.u1() > 0.0f ? RotB.col(1) : -RotB.col(1);
    frontNormal = -normal;
    front = posB * frontNormal + hB.u1();
    sideNormal = RotB.col(0);
    float side = posB * sideNormal;
    negSide = -side + hB.u0();
    posSide = side + hB.u0();
    negEdge = EDGE2;
    posEdge = EDGE4;
    ComputeIncidentEdge(incidentEdge, hA, posA, RotA, frontNormal);
  } break;
  }

  // clip other face with 5 box planes (1 face plane, 4 edge planes)

  ClipVertex clipPoints1[2];
  ClipVertex clipPoints2[2];
  int np;

  // Clip to box side 1
  np = ClipSegmentToLine(clipPoints1, incidentEdge, -sideNormal, negSide, negEdge);

  if (np < 2)
    return 0;

  // Clip to negative box side 1
  np = ClipSegmentToLine(clipPoints2, clipPoints1, sideNormal, posSide, posEdge);

  if (np < 2)
    return 0;

  // Now clipPoints2 contains the clipping points.
  // Due to roundoff, it is possible that clipping removes all points.
  //ClipVertex* cv = clipPoints2;
  ClipVertex* cv = incedge2;

  int numContacts = 0;
  for (int i = 0; i < 2; ++i) {
    float separation = frontNormal * cv[i].v - front;

    if (separation <= 0) {
      contacts[numContacts].separation = separation;
      contacts[numContacts].normal = normal;
      // slide contact point onto reference face (easy to cull)
      contacts[numContacts].position = cv[i].v - separation * frontNormal;
      contacts[numContacts].feature = cv[i].fp;
      if (axis == FACE_B_X || axis == FACE_B_Y)
        Flip(contacts[numContacts].feature);
      ++numContacts;
    }
  }

  return numContacts;
}