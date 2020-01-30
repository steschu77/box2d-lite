#include <math.h>

#include "VecMath.h"

// ----------------------------------------------------------------------------
static void xchg(unsigned k, const float a[2][2], float b[2][2]);

// ============================================================================
x3d::matrix2x2 x3d::getZeroMatrix2x2()
{
  return matrix2x2(0.0f);
}

// ============================================================================
x3d::matrix2x2 x3d::getIdentityMatrix2x2()
{
  return matrix2x2(1.0f);
}

// ============================================================================
x3d::matrix2x2::matrix2x2(float x)
{
  float c = cosf(x), s = sinf(x);
  m[0][0] = c;
  m[0][1] = -s;
  m[1][0] = s;
  m[1][1] = c;
}

// ----------------------------------------------------------------------------
x3d::matrix2x2::matrix2x2(float x00, float x01, float x10, float x11)
{
  m[0][0] = x00;
  m[0][1] = x01;
  m[1][0] = x10;
  m[1][1] = x11;
}

// ----------------------------------------------------------------------------
x3d::matrix2x2::matrix2x2(const float x[2][2])
{
  m[0][0] = x[0][0];
  m[0][1] = x[0][1];
  m[1][0] = x[1][0];
  m[1][1] = x[1][1];
}

// ----------------------------------------------------------------------------
x3d::matrix2x2 x3d::matrix2x2::transpose() const
{
  float a[2][2];

  for (int j = 0; j < 2; j++) {
    for (int i = 0; i < 2; i++) {
      a[j][i] = m[i][j];
    }
  }

  return x3d::matrix2x2(a);
}

#if 0
// ----------------------------------------------------------------------------
x3d::matrix2x2 x3d::matrix2x2::inverse() const
{
    float a = m[0][0];
    float b = m[0][1];
    float c = m[1][0];
    float d = m[1][1];

    float det = a * d - b * c;
    //assert(det != 0.0f);

    det = 1.0f / det;
    float x[2][2];
    x[0][0] = det * d;
    x[0][1] = -det * b;
    x[1][0] = -det * c;
    x[1][1] = det * a;

    return x3d::matrix2x2(x);
}
#endif
#if 1
// ----------------------------------------------------------------------------
x3d::matrix2x2 x3d::matrix2x2::inverse() const
{
  float a[2][2];
  float b[2][2];

  if (m[0][0] < epsilon || m[1][1] < epsilon) {
    return x3d::matrix2x2();
  }

  xchg(0, m, b);
  xchg(1, b, a);

  return x3d::matrix2x2(a);
}
#endif

// ----------------------------------------------------------------------------
static void xchg(unsigned k, const float a[2][2], float b[2][2])
{
  // (b)ij =      1 / (a)kk for i = j = k
  // (b)ij =  (a)ij / (a)kk for i = k, all j = 1..4, j != k
  // (b)ij = -(a)ij / (a)kk for j = k, all i = 1..4, i != k
  // (b)ij =  (a)ij - (a)ik * (a)kj / (a)kk for all i != k, j != k

  for (unsigned i = 0; i < 2; i++) {
    if (i == k) {
      for (unsigned j = 0; j < 2; j++) {
        if (j == k)
          b[i][j] = 1 / a[i][j];
        else
          b[i][j] = a[i][j] / a[k][k];
      }
    } else {
      for (unsigned j = 0; j < 2; j++) {
        if (j == k)
          b[i][j] = -a[i][j] / a[k][k];
        else
          b[i][j] = a[i][j] - a[i][k] * a[k][j] / a[k][k];
      }
    }
  }
}

// ----------------------------------------------------------------------------
// Solve M * x = v, where b is a column vector. This is more efficient
// than computing the inverse in one-shot cases.
x3d::vector2 x3d::matrix2x2::solve(const x3d::vector2& v) const
{
  float det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
  if (det != 0.0f) {
    det = 1.0f / det;
  }

  return x3d::vector2(det * (m[1][1] * v.m[0] - m[0][1] * v.m[1]),
    det * (m[0][0] * v.m[1] - m[1][0] * v.m[0]));
}

// ============================================================================
x3d::vector2 x3d::operator*(const x3d::vector2& v1, const x3d::matrix2x2& m2)
{
  float a0 = v1.m[0] * m2.m[0][0] + v1.m[1] * m2.m[1][0];
  float a1 = v1.m[0] * m2.m[0][1] + v1.m[1] * m2.m[1][1];

  return x3d::vector2(a0, a1);
}

// ============================================================================
x3d::vector2 x3d::operator*(const x3d::matrix2x2& m2, const x3d::vector2& v1)
{
  float a0 = m2.m[0][0] * v1.m[0] + m2.m[0][1] * v1.m[1];
  float a1 = m2.m[1][0] * v1.m[0] + m2.m[1][1] * v1.m[1];

  return x3d::vector2(a0, a1);
}

// ============================================================================
x3d::matrix2x2 x3d::operator*(const x3d::matrix2x2& m1, const x3d::matrix2x2& m2)
{
  float a00 = m1.m[0][0] * m2.m[0][0] + m1.m[0][1] * m2.m[1][0];
  float a01 = m1.m[0][0] * m2.m[0][1] + m1.m[0][1] * m2.m[1][1];

  float a10 = m1.m[1][0] * m2.m[0][0] + m1.m[1][1] * m2.m[1][0];
  float a11 = m1.m[1][0] * m2.m[0][1] + m1.m[1][1] * m2.m[1][1];

  return x3d::matrix2x2(a00, a01, a10, a11);
}

// ============================================================================
x3d::matrix2x2 x3d::operator+(const x3d::matrix2x2& m1, const x3d::matrix2x2& m2)
{
  float a00 = m1.m[0][0] + m2.m[0][0];
  float a01 = m1.m[0][1] + m2.m[0][1];

  float a10 = m1.m[1][0] + m2.m[1][0];
  float a11 = m1.m[1][1] + m2.m[1][1];

  return x3d::matrix2x2(a00, a01, a10, a11);
}

// ----------------------------------------------------------------------------
x3d::matrix2x2 x3d::matrix2x2::abs() const
{
  return x3d::matrix2x2(
    fabsf(m[0][0]), fabsf(m[0][1]), fabsf(m[1][0]), fabsf(m[1][1]));
}

// ----------------------------------------------------------------------------
x3d::vector2 x3d::matrix2x2::col(unsigned idx) const
{
  return x3d::vector2(m[0][idx], m[1][idx]);
}
