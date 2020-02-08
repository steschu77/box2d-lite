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
  ClipVertex() {}

  ContactPointId id;
  x3d::vector2 v;
};

struct ReferenceEdge
{
  ReferenceEdge(const Body* poly1, const Body* poly2, int flip);

  const Body* poly1 = nullptr;
  const Body* poly2 = nullptr;
  int flip = 0;
  float separation = -FLT_MAX;
  int index = 0;
};

ReferenceEdge::ReferenceEdge(const Body* poly1, const Body* poly2, int flip)
: poly1(poly1)
, poly2(poly2)
, flip(flip)
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
  c[0].id.setInEdge2(edge->index);
  c[0].id.setOutEdge2(i1);

  c[1].v = poly2->wVerts[i2];
  c[1].id.setInEdge2(edge->index);
  c[1].id.setOutEdge2(i2);
}

void ClipSegmentToLine(ClipVertex vOut[2],
  const x3d::vector2& normal, const x3d::vector2& vx, int clipEdge)
{
  // Calculate the distance of end points to the line
  float distance0 = normal * (vOut[0].v - vx) - 0.02f;
  float distance1 = normal * (vOut[1].v - vx) - 0.02f;

  if (distance0 > 0.0f) {
    float interp = distance0 / (distance0 - distance1);
    vOut[0].v = vOut[0].v + interp * (vOut[1].v - vOut[0].v);
    vOut[0].id.setInEdge1(clipEdge);
    vOut[0].id.resetInEdge2();
  } else if (distance1 > 0.0f) {
    float interp = distance0 / (distance0 - distance1);
    vOut[1].v = vOut[0].v + interp * (vOut[1].v - vOut[0].v);
    vOut[1].id.setOutEdge1(clipEdge);
    vOut[1].id.resetOutEdge2();
  }
}

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

  const x3d::vector2 normal = edge->poly1->wNorms[iv1];
  const x3d::vector2 tangent = -normal.perpendicular();

  ClipSegmentToLine(incedge, -tangent, v11, iv1);
  ClipSegmentToLine(incedge, tangent, v12, iv2);

  // Now clipPoints2 contains the clipping points.
  // Due to roundoff, it is possible that clipping removes all points.
  int numContacts = 0;
  for (int i = 0; i < 2; ++i) {
    float separation = normal * (incedge[i].v - v11);

    if (separation <= 0) {
      contacts[numContacts].separation = separation;
      contacts[numContacts].normal = edge->flip ? -normal : normal;
      contacts[numContacts].id = edge->flip ? -incedge[i].id : incedge[i].id;
      // slide contact point onto reference face (easy to cull)
      contacts[numContacts].position = incedge[i].v - separation * normal;
      ++numContacts;
    }
  }

  return numContacts;
}