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

#include "Arbiter.h"
#include "RigidBody.h"
#include "World.h"

#include <algorithm>

inline float clamp(float a, float low, float high)
{
  return std::max(low, std::min(a, high));
}

ArbiterKey::ArbiterKey(RigidBody* b1, RigidBody* b2)
{
  if (b1 < b2) {
    body1 = b1;
    body2 = b2;
  } else {
    body1 = b2;
    body2 = b1;
  }
}

Arbiter::Arbiter(ArbiterKey& key)
: key(key)
{
  friction = sqrtf(key.body1->getFriction() * key.body2->getFriction());
}

// ----------------------------------------------------------------------------
void Arbiter::updateContacts(Contact* newContacts, int numNewContacts)
{
  // Store accumulated impulses
  for (int i = 0; i < numNewContacts; ++i) {
    for (int j = 0; j < numContacts; ++j) {
      if (newContacts[i].id.value == contacts[j].id.value) {
        newContacts[i].Pn = contacts[j].Pn;
        newContacts[i].Pt = contacts[j].Pt;
        newContacts[i].Pnb = contacts[j].Pnb;
        break;
      }
    }
  }

  for (int i = 0; i < numNewContacts; ++i) {
    contacts[i] = newContacts[i];
  }

  numContacts = numNewContacts;
}

void Arbiter::PreStep(float inv_dt)
{
  RigidBody* b1 = key.body1;
  RigidBody* b2 = key.body2;

  const float k_allowedPenetration = 0.01f;
  const float k_biasFactor = 0.2f;

  for (int i = 0; i < numContacts; ++i) {
    Contact* c = contacts + i;

    x3d::vector2 r1 = c->position - b1->getPosition();
    x3d::vector2 r2 = c->position - b2->getPosition();

    // Precompute normal mass, tangent mass, and bias.
    float rn1 = r1 * c->normal;
    float rn2 = r2 * c->normal;
    float kNormal = b1->getInvMass() + b2->getInvMass();
    float kTangent = kNormal;
    kNormal += b1->getInvI() * (r1 * r1 - rn1 * rn1)
      + b2->getInvI() * (r2 * r2 - rn2 * rn2);
    c->massNormal = 1.0f / kNormal;

    x3d::vector2 tangent = -c->normal.perpendicular();
    float rt1 = r1 * tangent;
    float rt2 = r2 * tangent;
    kTangent += b1->getInvI() * (r1 * r1 - rt1 * rt1) + b2->getInvI() * (r2 * r2 - rt2 * rt2);
    c->massTangent = 1.0f / kTangent;

    c->bias = -k_biasFactor * inv_dt
      * std::min(0.0f, c->separation + k_allowedPenetration);

    // Apply normal + friction impulse
    x3d::vector2 P = c->Pn * c->normal + c->Pt * tangent;

    b1->applyImpulse(c->position, -P);
    b2->applyImpulse(c->position, P);
  }
}

void Arbiter::ApplyImpulse()
{
  RigidBody* b1 = key.body1;
  RigidBody* b2 = key.body2;

  for (int i = 0; i < numContacts; ++i) {
    Contact* c = contacts + i;

    // Relative velocity at contact
    x3d::vector2 dv = b2->relativeVelocity(c->position)
      - b1->relativeVelocity(c->position);

    // Compute normal impulse
    float vn = dv * c->normal;

    float dPn = c->massNormal * (-vn + c->bias);

    // Clamp the accumulated impulse
    float Pn0 = c->Pn;
    c->Pn = std::max(Pn0 + dPn, 0.0f);
    dPn = c->Pn - Pn0;

    // Relative velocity at contact
    dv = b2->relativeVelocity(c->position)
      - b1->relativeVelocity(c->position);

    x3d::vector2 tangent = -c->normal.perpendicular();
    float vt = dv * tangent;
    float dPt = c->massTangent * (-vt);

    // Compute friction impulse
    float maxPt = friction * c->Pn;

    // Clamp friction
    float oldTangentImpulse = c->Pt;
    c->Pt = clamp(oldTangentImpulse + dPt, -maxPt, maxPt);
    dPt = c->Pt - oldTangentImpulse;

    // Apply contact impulse
    x3d::vector2 Pn = dPn * c->normal;
    x3d::vector2 Pt = dPt * tangent;
    x3d::vector2 P = Pt + Pn;

    b1->applyImpulse(c->position, -P);
    b2->applyImpulse(c->position, P);
  }
}
