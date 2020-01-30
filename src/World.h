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

#ifndef WORLD_H
#define WORLD_H

#include "Arbiter.h"
#include "VecMath.h"
#include <map>
#include <set>
#include <vector>

struct Body;
struct Joint;

struct World
{
  World(x3d::vector2 gravity, int iterations)
  : gravity(gravity)
  , iterations(iterations)
  {
  }

  void Add(Body* body);
  void AddStatic(Body* body);
  void Add(Joint* joint);
  void Clear();

  void Step(float dt);

  std::vector<Body*> bodies;
  std::vector<Body*> statics;
  std::vector<Joint*> joints;
  std::map<ArbiterKey, Arbiter> arbiters;
  x3d::vector2 gravity;
  int iterations;

private:
  void BroadPhase();
  void _collide(Body* bi, Body* bj);
};

#endif
