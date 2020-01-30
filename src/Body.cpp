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

#include "Body.h"

void Body::Set(const x3d::vector2& pos, float rot, const x3d::vector2& w, float m)
{
  position = pos;
  rotation = rot;
  velocity = x3d::vector2(0.0f, 0.0f);
  angularVelocity = 0.0f;
  force = x3d::vector2(0.0f, 0.0f);
  torque = 0.0f;
  friction = 0.2f;

  width = w;
  mass = m;

  if (mass < FLT_MAX) {
    invMass = 1.0f / mass;
    I = mass * (width.u0() * width.u0() + width.u1() * width.u1()) / 12.0f;
    invI = 1.0f / I;
  } else {
    invMass = 0.0f;
    I = FLT_MAX;
    invI = 0.0f;
  }
}

void Body::SetStatic(const x3d::vector2& pos, float rot, const x3d::vector2& w)
{
  position = pos;
  rotation = rot;
  velocity = x3d::vector2(0.0f, 0.0f);
  angularVelocity = 0.0f;
  force = x3d::vector2(0.0f, 0.0f);
  torque = 0.0f;
  friction = 0.2f;

  width = w;
  mass = FLT_MAX;
  invMass = 0.0f;
  I = FLT_MAX;
  invI = 0.0f;
}

void Body::applyImpulse(const x3d::vector2& pt, const x3d::vector2& P)
{
  x3d::vector2 r = pt - position;
  velocity += invMass * P;
  angularVelocity += invI * x3d::cross(r, P);
}