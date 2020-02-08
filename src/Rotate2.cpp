#include <math.h>

#include "VecMath.h"

// ============================================================================
x3d::rot2 x3d::rotation2(float angle)
{
  float s = sinf(angle);
  float c = cosf(angle);
  return x3d::rot2(s, c);
}

// ============================================================================
x3d::rot2::rot2()
: s(0.0f)
, c(1.0f)
{
}

// ----------------------------------------------------------------------------
x3d::rot2::rot2(float s, float c)
: s(s)
, c(c)
{
}

// ----------------------------------------------------------------------------
x3d::rot2 x3d::mul(const x3d::rot2& q, const x3d::rot2& r)
{
  float s = q.s * r.c + q.c * r.s;
  float c = q.c * r.c - q.s * r.s;
  return x3d::rot2(s, c);
}

// ----------------------------------------------------------------------------
x3d::rot2 x3d::mulT(const x3d::rot2& q, const x3d::rot2& r)
{
  float s = q.c * r.s - q.s * r.c;
  float c = q.c * r.c + q.s * r.s;
  return x3d::rot2(s, c);
}

// Rotate a vector
// ----------------------------------------------------------------------------
x3d::vector2 x3d::mul(const x3d::rot2& q, const x3d::vector2& v)
{
  float u0 = q.c * v.u0() - q.s * v.u1();
  float u1 = q.s * v.u0() + q.c * v.u1();
  return x3d::vector2(u0, u1);
}

// Inverse rotate a vector
// ----------------------------------------------------------------------------
x3d::vector2 x3d::mulT(const x3d::rot2& q, const x3d::vector2& v)
{
  float u0 = q.c * v.u0() + q.s * v.u1();
  float u1 = -q.s * v.u0() + q.c * v.u1();
  return x3d::vector2(u0, u1);
}
