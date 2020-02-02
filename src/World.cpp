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

#include "World.h"
#include "Body.h"
#include "Joint.h"

typedef std::map<ArbiterKey, Arbiter>::iterator ArbIter;
typedef std::pair<ArbiterKey, Arbiter> ArbPair;

void World::Add(Body* body)
{
  bodies.push_back(body);
}

void World::AddStatic(Body* body)
{
  statics.push_back(body);
}

void World::Add(Joint* joint)
{
  joints.push_back(joint);
}

void World::Clear()
{
  bodies.clear();
  statics.clear();
  joints.clear();
  arbiters.clear();
}

void World::_collide(Body* bi, Body* bj)
{
  ArbiterKey key(bi, bj);

  Contact contacts[2];
  int numContacts = Collide(contacts, key.body1, key.body2);

  if (numContacts > 0) {
    std::pair<ArbIter, bool> ib = arbiters.insert(ArbPair(key, Arbiter(key)));
    ib.first->second.updateContacts(contacts, numContacts);
  } else {
    arbiters.erase(key);
  }
}

void World::BroadPhase()
{
  // O(n^2) broad-phase
  for (size_t i = 0; i < bodies.size(); ++i) {
    Body* bi = bodies[i];

    for (size_t j = i + 1; j < bodies.size(); ++j) {
      Body* bj = bodies[j];
      _collide(bi, bj);
    }

    for (size_t j = 0; j < statics.size(); ++j) {
      Body* bj = statics[j];
      _collide(bi, bj);
    }
  }
}

void World::Step(float dt)
{
  float inv_dt = dt > 0.0f ? 1.0f / dt : 0.0f;

  // Determine overlapping bodies and update contact points.
  BroadPhase();

  // Integrate forces.
  for (auto& b : bodies) {
    b->velocity += dt * (gravity + b->invMass * b->force);
    b->angularVelocity += dt * b->invI * b->torque;
  }

  // Perform pre-steps.
  for (auto& arb : arbiters) {
    arb.second.PreStep(inv_dt);
  }

  for (auto& joint : joints) {
    joint->PreStep(inv_dt);
  }

  // Perform iterations
  for (int i = 0; i < iterations; ++i) {
    for (auto& arb : arbiters) {
      arb.second.ApplyImpulse();
    }

    for (auto& joint : joints) {
      joint->ApplyImpulse();
    }
  }

  // Integrate Velocities
  for (auto& b : bodies) {
    b->position += dt * b->velocity;
    b->rotation += dt * b->angularVelocity;

    b->force = x3d::vector2(0.0f, 0.0f);
    b->torque = 0.0f;
  }
}
