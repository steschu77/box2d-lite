/*
 * Copyright (c) 2006-2009 Erin Catto http://www.gphysics.com
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies.
 * Erin Catto makes no representations about the suitability
 * of this software for any purpose.
 * It is provided "as is" without express or implied warranty.
 */

#ifndef ARBITER_H
#define ARBITER_H

#include "VecMath.h"

struct RigidBody;

struct ContactPointId
{
  ContactPointId() {}

  explicit ContactPointId(int val)
  : value(val)
  {
  }

  void setInEdge1(int val) { value |= val << 0; }
  void setOutEdge1(int val) { value |= val << 8; }
  void setInEdge2(int val) { value |= val << 16; }
  void setOutEdge2(int val) { value |= val << 24; }
  void resetInEdge2() { value &= ~0xff0000; }
  void resetOutEdge2() { value &= ~0xff000000; }

  ContactPointId operator-() const
  {
    return ContactPointId(((value & 0xffff) << 16) | ((value & 0xffff0000) >> 16));
  }

  int value = 0;
};

struct Contact
{
  x3d::vector2 position;
  x3d::vector2 normal;
  float separation = 0;
  float Pn = 0; // accumulated normal impulse
  float Pt = 0; // accumulated tangent impulse
  float Pnb = 0; // accumulated normal impulse for position bias
  float massNormal = 0;
  float massTangent = 0;
  float bias = 0;
  ContactPointId id;
};

struct ContactPoint
{
  ContactPointId id;
  x3d::vector2 v;
  x3d::vector2 normal;
  float separation = 0;
};

struct ContactPoints
{
  ContactPoint pt[2];
  int numContacts;
};

struct ArbiterKey
{
  ArbiterKey(RigidBody* b1, RigidBody* b2);

  RigidBody* body1;
  RigidBody* body2;
};

struct Arbiter
{
  enum
  {
    MAX_POINTS = 2
  };

  Arbiter(ArbiterKey& key);

  void updateContacts(Contact* contacts, int numContacts);

  void PreStep(float inv_dt);
  void ApplyImpulse();

  ArbiterKey key;

  Contact contacts[MAX_POINTS];
  int numContacts = 0;

  // Combined friction
  float friction = 0;
};

inline bool operator<(const ArbiterKey& x0, const ArbiterKey& x1)
{
  if (x0.body1 < x1.body1) {
    return true;
  }
  if (x0.body1 == x1.body1 && x0.body2 < x1.body2) {
    return true;
  }
  return false;
}

int Collide(Contact* contacts, RigidBody* body1, RigidBody* body2);

#endif
