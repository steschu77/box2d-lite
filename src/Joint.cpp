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
#include "Body.h"
#include "World.h"

void Joint::Set(Body* b1, Body* b2, const x3d::vector2& anchor)
{
  body1 = b1;
  body2 = b2;

  x3d::matrix2x2 Rot1(body1->rotation);
  x3d::matrix2x2 Rot2(body2->rotation);
  x3d::matrix2x2 Rot1T = Rot1.transpose();
  x3d::matrix2x2 Rot2T = Rot2.transpose();

  localAnchor1 = Rot1T * (anchor - body1->position);
  localAnchor2 = Rot2T * (anchor - body2->position);

  P = x3d::vector2(0.0f, 0.0f);

  softness = 0.0f;
  biasFactor = 0.2f;
}

void Joint::PreStep(float inv_dt)
{
  // Pre-compute anchors, mass matrix, and bias.
  x3d::matrix2x2 Rot1(body1->rotation);
  x3d::matrix2x2 Rot2(body2->rotation);

  r1 = Rot1 * localAnchor1;
  r2 = Rot2 * localAnchor2;

  // deltaV = deltaV0 + K * impulse
  // invM = [(1/m1 + 1/m2) * eye(2) - skew(r1) * invI1 * skew(r1) - skew(r2) *
  // invI2 * skew(r2)]
  //      = [1/m1+1/m2     0    ] + invI1 * [r1.u1()*r1.u1() -r1.u0()*r1.u1()] +
  //      invI2 * [r1.u1()*r1.u1() -r1.u0()*r1.u1()]
  //        [    0     1/m1+1/m2]           [-r1.u0()*r1.u1() r1.u0()*r1.u0()]
  //        [-r1.u0()*r1.u1() r1.u0()*r1.u0()]
  x3d::matrix2x2 K1;
  K1.m[0][0] = body1->invMass + body2->invMass;
  K1.m[0][1] = 0.0f;
  K1.m[1][0] = 0.0f;
  K1.m[1][1] = body1->invMass + body2->invMass;

  x3d::matrix2x2 K2;
  K2.m[0][0] = body1->invI * r1.u1() * r1.u1();
  K2.m[0][1] = -body1->invI * r1.u0() * r1.u1();
  K2.m[1][0] = -body1->invI * r1.u0() * r1.u1();
  K2.m[1][1] = body1->invI * r1.u0() * r1.u0();

  x3d::matrix2x2 K3;
  K3.m[0][0] = body2->invI * r2.u1() * r2.u1();
  K3.m[0][1] = -body2->invI * r2.u0() * r2.u1();
  K3.m[1][0] = -body2->invI * r2.u0() * r2.u1();
  K3.m[1][1] = body2->invI * r2.u0() * r2.u0();

  x3d::matrix2x2 K = K1 + K2 + K3;
  K.m[0][0] += softness;
  K.m[1][1] += softness;

  M = K.inverse();

  x3d::vector2 p1 = body1->position + r1;
  x3d::vector2 p2 = body2->position + r2;
  x3d::vector2 dp = p2 - p1;

  bias = -biasFactor * inv_dt * dp;

  // Apply accumulated impulse.
  body1->velocity -= body1->invMass * P;
  body1->angularVelocity -= body1->invI * x3d::cross(r1, P);

  body2->velocity += body2->invMass * P;
  body2->angularVelocity += body2->invI * x3d::cross(r2, P);
}

void Joint::ApplyImpulse()
{
  x3d::vector2 dv = body2->velocity + x3d::cross(body2->angularVelocity, r2)
    - body1->velocity - x3d::cross(body1->angularVelocity, r1);

  x3d::vector2 impulse;

  impulse = M * (bias - dv - softness * P);

  body1->velocity -= body1->invMass * impulse;
  body1->angularVelocity -= body1->invI * x3d::cross(r1, impulse);

  body2->velocity += body2->invMass * impulse;
  body2->angularVelocity += body2->invI * x3d::cross(r2, impulse);

  P += impulse;
}
