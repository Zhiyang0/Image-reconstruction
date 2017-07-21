// Andrea Gutierrez
// University College London
// Medical Physics and Biomedical Engineering Department
// November 2016

#ifndef GAMMAEVENT_H
#define GAMMAEVENT_H

#include <iostream>
#include <fstream>
#include <list>
#include <math.h>
#include <vector>
#include "VoxelParam.hh"
#include "Vec3.h"


#define X_POS 0
#define Y_POS 1
#define Z_POS 2

#define ELECTRON_MASS 510.999 // electron mass in keV
#define AXES 3

using namespace std;
class C3Vector;
class VoxelParam;
class FOVcalculator;
class GammaEvent{
public:
  GammaEvent(  float x1, float y1, float z1, float e1, float x2, float y2, float z2, float e2) : mscatt_x1(x1), mscatt_y1(y1), mscatt_z1(z1), mscatt_e1(e1), mabs_x2(x2), mabs_y2(y2), mabs_z2(z2), mabs_e2(e2),meventID(0),mdir(AXES), morg(AXES),rotZY(AXES,std::vector<float>(AXES,0)),newPosX(0),newPosY(0),newPosZ(0),oldPosX(0),oldPosY(0),oldPosZ(0),mPhiMin(0),mPhiMax(0) {}
  GammaEvent();
  ~GammaEvent(){};
  float GetScatterX(void) { return mscatt_x1; }
  float GetScatterY(void) { return mscatt_y1; }
  float GetScatterZ(void) { return mscatt_z1; }

  float GetAbsX(void) { return mabs_x2; }
  float GetAbsY(void) { return mabs_y2; }
  float GetAbsZ(void) { return mabs_z2; }
  
  float GetScatterE(void) {return mscatt_e1;}
  float GetAbsE(void) {return mabs_e2;}


  std::vector<float> GetDirectionVector(){ return mdir;}
  std::vector<float> GetOriginVector(){ return morg;}


  void SetPhiMax(float phi){
    mPhiMax = phi;
  }
  float GetPhiMax(){return mPhiMax;}

  void SetPhiMin(float phi){
    mPhiMin = phi;
  }
  float GetPhiMin(){return mPhiMin;}

  int GetEventID(void){return meventID;}
 
  void SetEventID(int event){
    meventID = event;
  }

  void SetNewPos(std::vector <float> &vec){
    newPosX = vec[X_POS];
    newPosY = vec[Y_POS];
    newPosZ = vec[Z_POS];
}

  void SetOldPos(std::vector <float> &vec){
    oldPosX = vec[X_POS];
    oldPosY = vec[Y_POS];
    oldPosZ = vec[Z_POS];
}
  float GetNewX(){return newPosX;}
  float GetNewY(){return newPosY;}
  float GetNewZ(){return newPosZ;}

  float GetOldX(){return oldPosX;}
  float GetOldY(){return oldPosY;}
  float GetOldZ(){return oldPosZ;}

  float GetRotY(){return rotAngY;}
  float GetRotZ(){return rotAngZ;}

  float GetComptonAngle(void) {return mcomptonAngle; }
  void CalculateComptonAngle();

  void CalculateConicalSurface();
  void CalculateNewRandomPosition(std::vector<float> &new_pos, VoxelParam * image);

  friend void VoxelParam::CalculateFOV(GammaEvent *event);

private:
  float mscatt_x1,mscatt_y1,mscatt_z1,mscatt_e1,mabs_x2,mabs_y2,mabs_z2,mabs_e2;
  float mcomptonAngle;
  int meventID;
  float newPosX, newPosY, newPosZ;
  float oldPosX, oldPosY, oldPosZ;
  float mPhiMin, mPhiMax;
 
  //cone Parameters
  std::vector<float> mdir; //axis direction
  std::vector<float> morg; //tip of the cone.
  float rotAngY,rotAngZ;
  //  std::vector<std::vector<float> > rotZY(AXES, std::vector<float>(AXES, 0));
  std::vector<std::vector<float> > rotZY;
  //(AXES, std::vector<float>(AXES, 0));
  
};

#endif
