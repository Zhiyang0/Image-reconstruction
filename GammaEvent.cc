// Andrea Gutierrez
// University College London
// Medical Physics and Biomedical Engineering Department
// November 2016

#include "GammaEvent.hh"
#include <stdlib.h>
#include "TRandom3.h"

#define X_POS 0
#define Y_POS 1
#define Z_POS 2
#define _USE_MATH_DEFINES
#include <cmath>
class Vec3;

float rand_FloatRange(float a, float b)
{
  return ((b-a)*((float)rand()/RAND_MAX))+a;
}

void GammaEvent::CalculateComptonAngle(){
  mcomptonAngle = acos(1 - ELECTRON_MASS*((1.0/mabs_e2) - 1.0/(mscatt_e1+mabs_e2)));
  //std::cout << "mcompton angle: " << mcomptonAngle << std::endl;
}


void GammaEvent::CalculateConicalSurface(){

  mdir[X_POS] = mscatt_x1 - mabs_x2;
  mdir[Y_POS] = mscatt_y1 - mabs_y2;
  mdir[Z_POS] = mscatt_z1 - mabs_z2;
  
  
  morg[X_POS] = mscatt_x1;
  morg[Y_POS] = mscatt_y1;
  morg[Z_POS] = mscatt_z1;
 

//float rotAngY, rotAngZ;

//calculating rotation Angle around Y    
rotAngY = acos((mdir[Z_POS]*1.0)/sqrt(mdir[X_POS]*mdir[X_POS] +mdir[Y_POS]*mdir[Y_POS] + mdir[Z_POS]*mdir[Z_POS]));

if (mdir[X_POS] < 0) 
  rotAngY = -1*rotAngY;


//rotation angle Z

rotAngZ = acos (mdir[X_POS]/sqrt(mdir[X_POS]*mdir[X_POS] + mdir[Y_POS]*mdir[Y_POS]));

 if((mdir[X_POS]>=0 && mdir[Y_POS]<0) || (mdir[X_POS]<0 && mdir[Y_POS]>=0))
  rotAngZ = -1*rotAngZ;
 
 //std::cout << "rot angle Y: " << rotAngY*180.0/PI << "\t" << rotAngZ*180.0/PI << std::endl;
 
//test: rotation matrix
 
 double angleMatY;
 double angleMatZ;

 Vec3 * axisZ = new Vec3( 0.0, 0.0, 1.0);
 Vec3 * axisX = new Vec3( 1.0, 0.0, 0.0);
 Vec3 * m_comptonAxisDirection = new Vec3( mdir[X_POS],mdir[Y_POS],mdir[Z_POS]);


 // Rotation around Y-axis, to rotate z_axis to X component of Compton axis
 angleMatY = m_comptonAxisDirection->GetScalarProductAngle(*axisZ);
 if (m_comptonAxisDirection->GetX() < 0)
   {
     angleMatY = -1.0 * angleMatY;
   }

 Vec3 * xyVector= new Vec3( fabs(m_comptonAxisDirection->GetX()), fabs(m_comptonAxisDirection->GetY()), 0.0);
 angleMatZ = xyVector->GetScalarProductAngle(*axisX);
 if (   (m_comptonAxisDirection->GetX() >= 0 && m_comptonAxisDirection->GetY() <  0)
	|| (m_comptonAxisDirection->GetX() <  0 && m_comptonAxisDirection->GetY() >= 0) )
   {
     angleMatZ = -1.0 * angleMatZ;
   }
 // std::cout << angleMatY*180.0/PI << "\t" <<angleMatZ*180.0/PI << std::endl; 

 rotAngY = angleMatY;
 rotAngZ = angleMatZ;


 
 //FillRotZYMatrix(rotZY,rotAngY,rotAngZ);
 rotZY[0][0] = cos(rotAngZ)*cos(rotAngY);
 rotZY[1][0] =-sin(rotAngZ);
 rotZY[2][0] = sin(rotAngY)*cos(rotAngZ);
 
 rotZY[0][1] = sin(rotAngZ)*cos(rotAngY);
 rotZY[1][1] = cos(rotAngZ);
 rotZY[2][1] = sin(rotAngZ)*sin(rotAngY);

 rotZY[0][2] = -sin(rotAngY);
 rotZY[1][2] = 0.0;
 rotZY[2][2] = cos(rotAngY);
 

}


void GammaEvent::CalculateNewRandomPosition(std::vector<float> &new_pos, VoxelParam * image){
  //get random Z
  float max_z = image->GetMaxZ();
  float min_z = image->GetMinZ();
 
  float xmin = image->GetMinX();// mm
  float xmax = image->GetMaxX();// mm
  float ymin = image->GetMinY();// mm
  float ymax = image->GetMaxY();// mm

  std::vector<float> cone(AXES);
  std::vector<float> rot_cone(AXES);
  float coeff, random_z, random_phi;
 
  //TRandom3 *r = new TRandom3(0);

  new_pos[X_POS]= -2000;
  new_pos[Y_POS]= -2000;
  int tries = 0;
   while ((tries < 1) && (new_pos[X_POS]<xmin || new_pos[X_POS]>xmax ||new_pos[Y_POS]<ymin || new_pos[Y_POS] > ymax)){
     if(image->GetZbins()==1){
       random_z = min_z;
     }
     else{

       random_z   = rand_FloatRange(min_z,max_z);//r->Uniform(min_z, max_z);//rand_FloatRange(min_z,max_z);

     }
     random_phi =rand_FloatRange(this->GetPhiMin(),this->GetPhiMax());//  r->Uniform(this->GetPhiMin(),this->GetPhiMax()); //rand_FloatRange(0,2*PI);
    
  cone[X_POS] = sin(mcomptonAngle)*cos(random_phi);
  cone[Y_POS] = sin(mcomptonAngle)*sin(random_phi); 
  cone[Z_POS] = cos(mcomptonAngle);
  
  rot_cone[X_POS]=rotZY[X_POS][X_POS]*cone[X_POS] + rotZY[Y_POS][X_POS]*cone[Y_POS]+ + rotZY[Z_POS][X_POS]*cone[Z_POS];
  rot_cone[Y_POS]=rotZY[X_POS][Y_POS]*cone[X_POS] + rotZY[Y_POS][Y_POS]*cone[Y_POS]+ + rotZY[Z_POS][Y_POS]*cone[Z_POS];
  rot_cone[Z_POS]=rotZY[X_POS][Z_POS]*cone[X_POS] + rotZY[Y_POS][Z_POS]*cone[Y_POS]+ + rotZY[Z_POS][Z_POS]*cone[Z_POS];

  coeff = (random_z-morg[Z_POS])/rot_cone[Z_POS];    
  
  new_pos[0] =  morg[X_POS]+coeff*rot_cone[X_POS];
  new_pos[1] =  morg[Y_POS]+coeff*rot_cone[Y_POS];
  new_pos[2] =  morg[Z_POS]+coeff*rot_cone[Z_POS];
  tries++;

   }
 
}


