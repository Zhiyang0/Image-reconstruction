#include "Vec3.h"


#include <cmath>


Vec3::Vec3(): m_x(0.0), m_y(0.0), m_z(0.0)
{
}

Vec3::Vec3(const double& in_x, const double& in_y, const double& in_z): m_x(in_x), m_y(in_y), m_z(in_z)
{
}

Vec3::~Vec3()
{
}


void
Vec3::Set (const double& in_x, const double& in_y, const double& in_z)
{
  m_x = in_x;
  m_y = in_y;
  m_z = in_z;
}

double
Vec3::GetLength() const
{
  double len = sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
  return len;
}

double
Vec3::GetScalarProduct(const Vec3& in_vec) const
{
  double inp = (*this) * in_vec;
  double len1 = GetLength();
  double len2 = in_vec.GetLength();
  
  inp = ((len1 * len2) != 0) ? (inp / (len1 * len2)) : 0;
  return inp;
}

double
Vec3::GetScalarProductAngle(const Vec3& in_vec) const
{
  double cosangle = GetScalarProduct(in_vec);
  
  double angle = acos(cosangle);
  return angle;
}


// + (plus) operator
Vec3 Vec3::operator+ (const Vec3& in_v) const
{
      Vec3 result;
      result.m_x = (m_x + in_v.GetX());
      result.m_y = (m_y + in_v.GetY());
      result.m_z = (m_z + in_v.GetZ());
      return result;
}

// - (minus) operator
Vec3 Vec3::operator- (const Vec3& in_v) const
{
      Vec3 result;
      result.m_x = (m_x - in_v.GetX());
      result.m_y = (m_y - in_v.GetY());
      result.m_z = (m_z - in_v.GetZ());
      return result;
}

// inproduct
double Vec3::operator*(const Vec3& in_v) const
{
	double result = m_x * in_v.GetX() + m_y * in_v.GetY() + m_z * in_v.GetZ();
	return result;
}

Vec3 Vec3::operator*(double const& in_coeff) const
{
	Vec3 result;
	result.m_x = (m_x * in_coeff);
	result.m_y = (m_y * in_coeff);
	result.m_z = (m_z * in_coeff);
	return result;
}




