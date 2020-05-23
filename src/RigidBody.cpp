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

#include "RigidBody.h"

// ----------------------------------------------------------------------------
RigidBody::RigidBody()
: _Position(0, 0)
, _LinearVelocity(0, 0)
, _Force(0, 0)
, _Rotation(0)
, _AngularVelocity(0)
, _Torque(0)
, _Mass(0)
, _invMass(0)
, _invInertia(0)
, _Friction(0.0f)
, _Restitution(0.0f)
, _Size(0, 0)
{
}

  // ----------------------------------------------------------------------------
RigidBody::RigidBody(int id, const x3d::vector2& pos, float rot, const x3d::vector2& w, float mass)
: _Position(pos)
, _LinearVelocity(0, 0)
, _Force(0, 0)
, _Rotation(rot)
, _AngularVelocity(0)
, _Torque(0)
, _Mass(mass)
, _invMass(0)
, _invInertia(0)
, _Friction(0.2f)
, _Restitution(0.0f)
, _Size(w)
{
  x3d::vector2 h = 0.5f * w;
  vertices[0] = x3d::vector2(-h.u0(), -h.u1());
  vertices[1] = x3d::vector2(h.u0(), -h.u1());
  vertices[2] = x3d::vector2(h.u0(), h.u1());
  vertices[3] = x3d::vector2(-h.u0(), h.u1());

  normals[0] = x3d::vector2(0.0f, -1.0f);
  normals[1] = x3d::vector2(1.0f, 0.0f);
  normals[2] = x3d::vector2(0.0f, 1.0f);
  normals[3] = x3d::vector2(-1.0f, 0.0f);

  if (mass < FLT_MAX) {
    _invMass = 1.0f / mass;
    const float I = mass * (w.u0() * w.u0() + w.u1() * w.u1()) / 12.0f;
    _invInertia = 1.0f / I;
  }

  updateWorld();
}

// ----------------------------------------------------------------------------
void RigidBody::Set(int id, const x3d::vector2& pos, float rot,
  const x3d::vector2& w, float mass)
{
  _Position = pos;
  _LinearVelocity = x3d::vector2(0, 0);
  _Force = x3d::vector2(0, 0);
  _Rotation = rot;
  _AngularVelocity = 0;
  _Torque = 0;
  _Mass = mass;
  _invMass = 0;
  _invInertia = 0;
  _Friction = 0.2f;
  _Restitution = 0.0f;
  _Size = w;

  x3d::vector2 h = 0.5f * w;
  vertices[0] = x3d::vector2(-h.u0(), -h.u1());
  vertices[1] = x3d::vector2(h.u0(), -h.u1());
  vertices[2] = x3d::vector2(h.u0(), h.u1());
  vertices[3] = x3d::vector2(-h.u0(), h.u1());

  normals[0] = x3d::vector2(0.0f, -1.0f);
  normals[1] = x3d::vector2(1.0f, 0.0f);
  normals[2] = x3d::vector2(0.0f, 1.0f);
  normals[3] = x3d::vector2(-1.0f, 0.0f);

  if (mass > 0) {
    const float I = mass * (w.u0() * w.u0() + w.u1() * w.u1()) / 12.0f;
    _invMass = 1.0f / mass;
    _invInertia = 1.0f / I;
  }

  updateWorld();
}

// ----------------------------------------------------------------------------
x3d::vector2 RigidBody::localPosition(const x3d::vector2& worldPt) const
{
  const x3d::vector2 r = worldPt - _Position;
  return r.rotate(-_Rotation);
}

// ----------------------------------------------------------------------------
x3d::vector2 RigidBody::worldPosition(const x3d::vector2& localPt) const
{
  return localPt.rotate(_Rotation) + _Position;
}

// ----------------------------------------------------------------------------
x3d::vector2 RigidBody::relativeVelocity(const x3d::vector2& pt) const
{
  return _LinearVelocity + x3d::cross(_AngularVelocity, pt - _Position);
}

// ----------------------------------------------------------------------------
void RigidBody::applyForce(const x3d::vector2& F)
{
  _Force += F;
}

// ----------------------------------------------------------------------------
void RigidBody::applyImpulse(const x3d::vector2& pt, const x3d::vector2& P)
{
  const x3d::vector2 r = pt - _Position;
  _LinearVelocity += _invMass * P;
  _AngularVelocity += _invInertia * x3d::cross(r, P);
}

// ----------------------------------------------------------------------------
void RigidBody::integrateForces(float dt)
{
  _LinearVelocity += dt * _invMass * _Force;
  _AngularVelocity += dt * _invInertia * _Torque;

  _Force = x3d::vector2(0.0f, 0.0f);
  _Torque = 0.0f;
}

// ----------------------------------------------------------------------------
void RigidBody::integrateVelocities(float dt)
{
  _Position += dt * _LinearVelocity;
  _Rotation += dt * _AngularVelocity;
  updateWorld();
}

// ----------------------------------------------------------------------------
void RigidBody::updateWorld()
{
  p = _Position;
  q = x3d::rotation2(_Rotation);

  for (int i = 0; i < count; i++) {
    wVerts[i] = x3d::mul(q, vertices[i]) + p;
  }

  for (int i = 0; i < count; i++) {
    wNorms[i] = x3d::mul(q, normals[i]);
  }
}
