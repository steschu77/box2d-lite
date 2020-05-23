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
#include "RigidBody.h"
#include "Joint.h"

typedef std::map<ArbiterKey, Arbiter>::iterator ArbIter;
typedef std::pair<ArbiterKey, Arbiter> ArbPair;

void World::Add(RigidBody* body)
{
  bodies.push_back(body);
}

void World::AddStatic(RigidBody* body)
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

void World::_collide(RigidBody* bi, RigidBody* bj)
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
    RigidBody* bi = bodies[i];

    for (size_t j = i + 1; j < bodies.size(); ++j) {
      RigidBody* bj = bodies[j];
      _collide(bi, bj);
    }

    for (size_t j = 0; j < statics.size(); ++j) {
      RigidBody* bj = statics[j];
      _collide(bi, bj);
    }
  }
}

void World::Step(float dt)
{
  float inv_dt = dt > 0.0f ? 1.0f / dt : 0.0f;

  // Determine overlapping bodies and update contact points.
  BroadPhase();

  for (auto& b : bodies) {
    b->applyForce(gravity * b->getMass());
  }

  // Integrate forces.
  for (auto& b : bodies) {
    b->integrateForces(dt);
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
    b->integrateVelocities(dt);
  }
}
