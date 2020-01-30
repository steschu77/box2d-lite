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

struct Body
{
  void Set(const x3d::vector2& pos, float rot, const x3d::vector2& w, float m);
  void SetStatic(const x3d::vector2& pos, float rot, const x3d::vector2& w);
  void applyImpulse(const x3d::vector2& pt, const x3d::vector2& P);

  x3d::vector2 position;
  float rotation = 0;

  x3d::vector2 velocity;
  float angularVelocity = 0;

  x3d::vector2 force;
  float torque = 0;

  x3d::vector2 width;

  float friction = 0.2f;
  float mass = FLT_MAX;
  float invMass = 0;
  float I = FLT_MAX;
  float invI = 0;
};

#endif
