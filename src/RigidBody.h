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

#ifndef BODY_H
#define BODY_H

#include "VecMath.h"

// ============================================================================
struct MassData
{
  MassData();
  MassData(float Mass, float Inertia);

  float Mass;
  float invMass;

  float Inertia;
  float invInertia;
};

// ============================================================================
struct MaterialData
{
  float density;
  float restitution;
  float staticFriction;
  float dynamicFriction;
};

// ============================================================================
struct RigidBody
{
  RigidBody();
  RigidBody(int id, const x3d::vector2& pos, float rot, const x3d::vector2& w, float mass);

  void Set(int id, const x3d::vector2& pos, float rot, const x3d::vector2& w, float mass);

  int getId() const;

  x3d::vector2 getPosition() const;
  x3d::vector2 getSize() const;
  float getRotation() const;

  float getFriction() const;
  void setFriction(float friction);

  float getMass() const;
  float getInvMass() const;
  float getInvI() const;

  x3d::vector2 getVelocity() const;
  float getAngularVelocity() const;

  x3d::vector2 localPosition(const x3d::vector2& worldPt) const;
  x3d::vector2 worldPosition(const x3d::vector2& localPt) const;
  x3d::vector2 relativeVelocity(const x3d::vector2& pt) const;

  void applyForce(const x3d::vector2& F);
  void applyImpulse(const x3d::vector2& pt, const x3d::vector2& P);
  void integrateForces(float dt);
  void integrateVelocities(float dt);

private:
  int _Id = 0;

  x3d::vector2 _Size;

  float _Mass;
  float _invMass;
  float _invInertia;

  float _Friction;
  float _Restitution;

  // Linear components
  x3d::vector2 _Position;
  x3d::vector2 _LinearVelocity;
  x3d::vector2 _Force;

  // Angular components
  float _Rotation; // radians
  float _AngularVelocity;
  float _Torque;

public:
  x3d::vector2 vertices[4];
  x3d::vector2 normals[4];
  const int count = 4;

  x3d::vector2 p;
  x3d::rot2 q;

  x3d::vector2 wVerts[4];
  x3d::vector2 wNorms[4];
  void updateWorld();
};

#endif

// ----------------------------------------------------------------------------
inline int RigidBody::getId() const
{
  return _Id;
}

// ----------------------------------------------------------------------------
inline x3d::vector2 RigidBody::getPosition() const
{
  return _Position;
}

// ----------------------------------------------------------------------------
inline x3d::vector2 RigidBody::getVelocity() const
{
  return _LinearVelocity;
}

// ----------------------------------------------------------------------------
inline float RigidBody::getAngularVelocity() const
{
  return _AngularVelocity;
}

// ----------------------------------------------------------------------------
inline x3d::vector2 RigidBody::getSize() const
{
  return _Size;
}

// ----------------------------------------------------------------------------
inline float RigidBody::getRotation() const
{
  return _Rotation;
}

// ----------------------------------------------------------------------------
inline float RigidBody::getFriction() const
{
  return _Friction;
}

// ----------------------------------------------------------------------------
inline void RigidBody::setFriction(float friction)
{
  _Friction = friction;
}

// ----------------------------------------------------------------------------
inline float RigidBody::getMass() const
{
  return _Mass;
}
  
// ----------------------------------------------------------------------------
inline float RigidBody::getInvMass() const
{
  return _invMass;
}

// ----------------------------------------------------------------------------
inline float RigidBody::getInvI() const
{
  return _invInertia;
}
