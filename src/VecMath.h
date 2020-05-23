#pragma once

#include <float.h>

namespace x3d {

extern const float epsilon;
extern const float k_pi;

// ----------------------------------------------------------------------------
struct vector2
{
  vector2();
  vector2(const vector2& v);

  vector2(float u0, float u1);
  vector2(const float u[2]);

  void get(float a[2]) const;
  void put(const float a[2]);

  float u0() const;
  float u1() const;

  vector2 operator=(const vector2& v);
  vector2 operator-() const;
  vector2 operator+=(const vector2& v);
  vector2 operator-=(const vector2& v);
  vector2 operator+(const vector2& v) const;
  vector2 operator-(const vector2& v) const;

  vector2 operator*(float f) const;
  vector2 operator/(float f) const;

  float operator*(const vector2& v) const;

  bool operator==(const vector2& v) const;
  bool operator!=(const vector2& v) const;

  float length2() const;
  float length() const;

  vector2 norm() const;
  vector2 perpendicular() const;
  vector2 rotate(float rad) const;

  vector2 abs() const;

  float m[2];
};

// ============================================================================
vector2 operator*(float f, const vector2& v);

// ============================================================================
struct vector3
// ----------------------------------------------------------------------------
{
  vector3(float u0 = 0, float u1 = 0, float u2 = 0);
  vector3(const float u[3]);
  vector3(const vector3& v);
  ~vector3();

  void get(float a[4]) const;
  void put(const float a[4]);

  float u0() const;
  float u1() const;
  float u2() const;

  vector3 operator=(const vector3& v);
  vector3 operator-() const;
  vector3 operator+=(const vector3& v);
  vector3 operator-=(const vector3& v);
  vector3 operator+(const vector3& v) const;
  vector3 operator-(const vector3& v) const;

  vector3 operator*(float f) const;
  vector3 operator/(float f) const;

  float operator*(const vector3& v) const;

  bool operator==(const vector3& v) const;
  bool operator!=(const vector3& v) const;

  float length() const;
  vector3 norm();
  vector3 cross(const vector3& v) const;

  float m[3];
};

vector3 operator*(float f, const vector3& v);

// ============================================================================
struct matrix2x2
// ----------------------------------------------------------------------------
{
  matrix2x2(float x = 0);
  matrix2x2(float x00, float x01, float x10, float x11);
  matrix2x2(const float x[2][2]);

  vector2 col(unsigned idx) const;

  matrix2x2 abs() const;

  matrix2x2 transpose() const;
  matrix2x2 inverse() const;
  vector2 solve(const vector2& v) const;

  float m[2][2];
};

vector2 operator*(const vector2& v, const matrix2x2& m);
vector2 operator*(const matrix2x2& m, const vector2& v);
matrix2x2 operator+(const matrix2x2& m1, const matrix2x2& m2);
matrix2x2 operator*(const matrix2x2& m1, const matrix2x2& m2);

// ============================================================================
struct matrix3x3
// ----------------------------------------------------------------------------
{
  matrix3x3(float x = 0);
  matrix3x3(float x00, float x01, float x02, float x10, float x11, float x12,
    float x20, float x21, float x22);
  matrix3x3(const float x[3][3]);
  ~matrix3x3();

  matrix3x3 transpose() const;
  matrix3x3 inverse() const;

  float m[3][3];
};

// ============================================================================
matrix2x2 getZeroMatrix2x2();
matrix2x2 getIdentityMatrix2x2();

matrix3x3 getZeroMatrix3x3();
matrix3x3 getIdentityMatrix3x3();

// ============================================================================
vector2 rotate2(const vector2& v, float deg);
matrix2x2 translate2x2(const vector2& v);
matrix2x2 rotate2x2(float deg);

float dot(const vector2& v0, const vector2& v1);
float dot(const vector3& v0, const vector3& v1);

float cross(const vector2& v0, const vector2& v1);
vector2 cross(const vector2& v, float s);
vector2 cross(float s, const vector2& v);
vector3 cross(const vector3& v0, const vector3& v1);
float distance(const vector2& x0, const vector2& x1);

// ----------------------------------------------------------------------------
// line3
x3d::vector3 line_from_points(const vector2& x0, const vector2& x1);
x3d::vector3 line_from_points(const vector3& x0, const vector3& x1);
x3d::vector3 intersect(const vector3& v0, const vector3& v1);
float gap(const vector3& line, const vector3& pt);

// ----------------------------------------------------------------------------
struct rot2
{
  rot2();
  rot2(float s, float c);

  /// Sine and cosine
  float s;
  float c;
};

rot2 rotation2(float Angle);
rot2 mul(const rot2& q, const rot2& r);
rot2 mulT(const rot2& q, const rot2& r);
vector2 mul(const rot2& q, const vector2& v);
vector2 mulT(const rot2& q, const vector2& v);

} // namespace x3d
