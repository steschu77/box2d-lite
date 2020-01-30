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

// The normal points from A to B
int Collide(Contact* contacts, Body* bodyA, Body* bodyB)
{
  // Setup
  x3d::vector2 hA = 0.5f * bodyA->width;
  x3d::vector2 hB = 0.5f * bodyB->width;

  x3d::vector2 posA = bodyA->position;
  x3d::vector2 posB = bodyB->position;

  x3d::matrix2x2 RotA(bodyA->rotation);
  x3d::matrix2x2 RotB(bodyB->rotation);

  x3d::matrix2x2 RotAT = RotA.transpose();
  x3d::matrix2x2 RotBT = RotB.transpose();

  x3d::vector2 a1 = RotA.col(0);
  x3d::vector2 a2 = RotA.col(1);
  x3d::vector2 b1 = RotB.col(0);
  x3d::vector2 b2 = RotB.col(1);

  x3d::vector2 dp = posB - posA;
  x3d::vector2 dA = RotAT * dp;
  x3d::vector2 dB = RotBT * dp;

  x3d::matrix2x2 C = RotAT * RotB;
  x3d::matrix2x2 absC = C.abs();
  x3d::matrix2x2 absCT = absC.transpose();

  // Box A faces
  x3d::vector2 faceA = dA.abs() - 1.01f * hA - absC * hB;
  x3d::vector2 faceB = dB.abs() - absCT * hA - 1.01f * hB;
  if (faceA.u0() > 0.0f || faceA.u1() > 0.0f || faceB.u0() > 0.0f
    || faceB.u1() > 0.0f) {
    return 0;
  }

  // Find best axis
  Axis axis;
  float separation;
  x3d::vector2 normal;

  // Box A faces
  axis = FACE_A_X;
  separation = faceA.u0();
  normal = dA.u0() > 0.0f ? RotA.col(0) : -RotA.col(0);

  const float relativeTol = 0.95f;
  const float absoluteTol = 0.01f;

  if (faceA.u1() > relativeTol * separation) {
    axis = FACE_A_Y;
    separation = faceA.u1();
    normal = dA.u1() > 0.0f ? RotA.col(1) : -RotA.col(1);
  }

  // Box B faces
  if (faceB.u0() > relativeTol * separation) {
    axis = FACE_B_X;
    separation = faceB.u0();
    normal = dB.u0() > 0.0f ? RotB.col(0) : -RotB.col(0);
  }

  if (faceB.u1() > relativeTol * separation) {
    axis = FACE_B_Y;
    separation = faceB.u1();
    normal = dB.u1() > 0.0f ? RotB.col(1) : -RotB.col(1);
  }

  // Setup clipping plane data based on the separating axis
  x3d::vector2 frontNormal, sideNormal;
  ClipVertex incidentEdge[2];
  float front, negSide, posSide;
  char negEdge, posEdge;

  // Compute the clipping lines and the line segment to be clipped.
  switch (axis) {
  case FACE_A_X: {
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

  int numContacts = 0;
  for (int i = 0; i < 2; ++i) {
    float separation = frontNormal * clipPoints2[i].v - front;

    if (separation <= 0) {
      contacts[numContacts].separation = separation;
      contacts[numContacts].normal = normal;
      // slide contact point onto reference face (easy to cull)
      contacts[numContacts].position = clipPoints2[i].v - separation * frontNormal;
      contacts[numContacts].feature = clipPoints2[i].fp;
      if (axis == FACE_B_X || axis == FACE_B_Y)
        Flip(contacts[numContacts].feature);
      ++numContacts;
    }
  }

  return numContacts;
}