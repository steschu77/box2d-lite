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
  const x3d::vector2* n1s = poly1->wNorms;
  const x3d::vector2* v1s = poly1->wVerts;
  const x3d::vector2* v2s = poly2->wVerts;

  for (int i = 0; i < count1; ++i) {

    // Get poly1 normal and vertex in poly2 coords.
    const x3d::vector2 n = n1s[i];
    const x3d::vector2 v1 = v1s[i];

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

  const x3d::vector2& n1 = poly1->wNorms[edge->index];
  const x3d::vector2* n2s = poly2->wNorms;

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

  c[0].v = poly2->wVerts[i1];
  c[0].fp.e.inEdge2 = edge->index;
  c[0].fp.e.outEdge2 = i1;

  c[1].v = poly2->wVerts[i2];
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
      vOut[numOut].fp.e.inEdge2 = 0;
    } else {
      vOut[numOut].fp = vIn[1].fp;
      vOut[numOut].fp.e.outEdge1 = clipEdge;
      vOut[numOut].fp.e.outEdge2 = 0;
    }
    ++numOut;
  }

  return numOut;
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
  const x3d::vector2* v1s = edge->poly1->wVerts;

  const int iv1 = edge->index;
  const int iv2 = iv1 + 1 < count1 ? iv1 + 1 : 0;

  const x3d::vector2 v11 = v1s[iv1];
  const x3d::vector2 v12 = v1s[iv2];

	const x3d::vector2 tangent = (v12 - v11).norm();
  const x3d::vector2 xnormal = tangent.perpendicular();
  const float frontOffset = xnormal * v11;
  const float sideOffset1 = -tangent * v11 + 0.02f;
  const float sideOffset2 = tangent * v12 + 0.02f;

	ClipVertex incedge1[2];
  int np1 = ClipSegmentToLine(incedge1, incedge, -tangent, sideOffset1, iv1);
  if (np1 < 2) {
    return 0;
  }

	ClipVertex incedge2[2];
  int np2 = ClipSegmentToLine(incedge2, incedge1, tangent, sideOffset2, iv2);
  if (np2 < 2) {
    return 0;
  }

  // Now clipPoints2 contains the clipping points.
  // Due to roundoff, it is possible that clipping removes all points.
  ClipVertex* cv = incedge2;

  int numContacts = 0;
  for (int i = 0; i < 2; ++i) {
    float separation = xnormal * cv[i].v - frontOffset;

    if (separation <= 0) {
      contacts[numContacts].separation = separation;
      contacts[numContacts].normal = edge->flip ? -xnormal : xnormal;
      // slide contact point onto reference face (easy to cull)
      contacts[numContacts].position = cv[i].v - separation * xnormal;
      contacts[numContacts].feature = cv[i].fp;
      if (edge->flip) {
        Flip(contacts[numContacts].feature);
      }
      ++numContacts;
    }
  }

  return numContacts;
}