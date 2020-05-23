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

#include "Joint.h"
#include "RigidBody.h"
#include "World.h"

void Joint::Set(RigidBody* b1, RigidBody* b2, const x3d::vector2& anchor)
{
  body1 = b1;
  body2 = b2;

  x3d::matrix2x2 Rot1(body1->getRotation());
  x3d::matrix2x2 Rot2(body2->getRotation());
  x3d::matrix2x2 Rot1T = Rot1.transpose();
  x3d::matrix2x2 Rot2T = Rot2.transpose();

  localAnchor1 = Rot1T * (anchor - body1->getPosition());
  localAnchor2 = Rot2T * (anchor - body2->getPosition());

  P = x3d::vector2(0.0f, 0.0f);

  softness = 0.0f;
  biasFactor = 0.2f;
}

void Joint::PreStep(float inv_dt)
{
  // Pre-compute anchors, mass matrix, and bias.
  x3d::matrix2x2 Rot1(body1->getRotation());
  x3d::matrix2x2 Rot2(body2->getRotation());

  r1 = Rot1 * localAnchor1;
  r2 = Rot2 * localAnchor2;

  float invMass1 = body1->getInvMass();
  float invMass2 = body2->getInvMass();

  float invI1 = body1->getInvI();
  float invI2 = body2->getInvI();

  // deltaV = deltaV0 + K * impulse
  // invM = [(1/m1 + 1/m2) * eye(2) - skew(r1) * invI1 * skew(r1) - skew(r2) *
  // invI2 * skew(r2)]
  //      = [1/m1+1/m2     0    ] + invI1 * [r1.u1()*r1.u1() -r1.u0()*r1.u1()] +
  //      invI2 * [r1.u1()*r1.u1() -r1.u0()*r1.u1()]
  //        [    0     1/m1+1/m2]           [-r1.u0()*r1.u1() r1.u0()*r1.u0()]
  //        [-r1.u0()*r1.u1() r1.u0()*r1.u0()]
  x3d::matrix2x2 K1;
  K1.m[0][0] = invMass1 + invMass2;
  K1.m[0][1] = 0.0f;
  K1.m[1][0] = 0.0f;
  K1.m[1][1] = invMass1 + invMass2;

  x3d::matrix2x2 K2;
  K2.m[0][0] = invI1 * r1.u1() * r1.u1();
  K2.m[0][1] = -invI1 * r1.u0() * r1.u1();
  K2.m[1][0] = -invI1 * r1.u0() * r1.u1();
  K2.m[1][1] = invI1 * r1.u0() * r1.u0();

  x3d::matrix2x2 K3;
  K3.m[0][0] = invI2 * r2.u1() * r2.u1();
  K3.m[0][1] = -invI2 * r2.u0() * r2.u1();
  K3.m[1][0] = -invI2 * r2.u0() * r2.u1();
  K3.m[1][1] = invI2 * r2.u0() * r2.u0();

  x3d::matrix2x2 K = K1 + K2 + K3;
  K.m[0][0] += softness;
  K.m[1][1] += softness;

  M = K.inverse();

  x3d::vector2 p1 = body1->getPosition() + r1;
  x3d::vector2 p2 = body2->getPosition() + r2;
  x3d::vector2 dp = p2 - p1;

  bias = -biasFactor * inv_dt * dp;

  // Apply accumulated impulse.
  body1->applyImpulse(p1, -P);
  body2->applyImpulse(p2, P);
}

void Joint::ApplyImpulse()
{
  x3d::vector2 p1 = body1->getPosition() + r1;
  x3d::vector2 p2 = body2->getPosition() + r2;

  x3d::vector2 dv = body2->relativeVelocity(p2) - body1->relativeVelocity(p1);

  x3d::vector2 impulse;

  impulse = M * (bias - dv - softness * P);
  body1->applyImpulse(p1, -impulse);
  body2->applyImpulse(p2, impulse);

  P += impulse;
}
