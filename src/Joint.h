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

#ifndef JOINT_H
#define JOINT_H

#include "VecMath.h"

struct RigidBody;

struct Joint
{
	Joint() :
		body1(0), body2(0),
		P(0.0f, 0.0f),
		biasFactor(0.2f), softness(0.0f)
		{}

	void Set(RigidBody* body1, RigidBody* body2, const x3d::vector2& anchor);

	void PreStep(float inv_dt);
	void ApplyImpulse();

	x3d::matrix2x2 M;
	x3d::vector2 localAnchor1, localAnchor2;
	x3d::vector2 r1, r2;
	x3d::vector2 bias;
	x3d::vector2 P;		// accumulated impulse
	RigidBody* body1;
	RigidBody* body2;
	float biasFactor;
	float softness;
};

#endif