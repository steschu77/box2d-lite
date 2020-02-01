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

struct Body;

union FeaturePair
{
	struct Edges
	{
		char inEdge1;
		char outEdge1;
		char inEdge2;
		char outEdge2;
	} e;
	int value;
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
	FeaturePair feature;
};

struct ArbiterKey
{
  ArbiterKey(Body* b1, Body* b2);

  Body* body1;
  Body* body2;
};

struct Arbiter
{
	enum {MAX_POINTS = 2};

	Arbiter(ArbiterKey& key);

	void Update(Contact* contacts, int numContacts);

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

int Collide(Contact* contacts, Body* body1, Body* body2);

#endif
