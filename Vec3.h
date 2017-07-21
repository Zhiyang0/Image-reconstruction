
#ifndef VEC3_H
#define VEC3_H

#include <iostream>
#include <fstream>

using namespace std;

class Vec3
{
public:
  Vec3();
  Vec3(const double& in_x, const double& in_y, const double& in_z);
  virtual ~Vec3();
  
 
  void Set (const double& in_x, const double& in_y, const double& in_z);
  inline double	GetX() const { return m_x; }
  inline double	GetY() const { return m_y; }
  inline double	GetZ() const { return m_z; }
  inline void SetX(const double& in_val) { m_x = in_val; }
  inline void SetY(const double& in_val) { m_y = in_val; }
  inline void SetZ(const double& in_val) { m_z = in_val; }
  

  double GetLength() const;
  double GetScalarProduct(const Vec3& in_vec) const;       // scalar product
  double GetScalarProductAngle(const Vec3& in_vec) const;  // scalar product angle
  

  Vec3 operator+(const Vec3&) const;       		// operator+()
  Vec3 operator-(const Vec3&) const;       		// operator-()
  
  double operator*(const Vec3& in_vector) const;     // inproduct
  Vec3 operator*(const double& in_coeff) const;      // multiplication factor
  
 private:
  // data
  double m_x, m_y, m_z;
	

};

#endif




