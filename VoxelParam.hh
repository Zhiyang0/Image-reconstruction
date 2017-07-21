// Andrea Gutierrez
// University College London
// Medical Physics and Biomedical Engineering Department
// November 2016

#ifndef VOXELPARAM_H
#define VOXELPARAM_H
class GammaEvent;
class IntersectionParam;
class VoxelParam
{

 public:
  VoxelParam(int Xbins, int Ybins, int Zbins, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, float dimx, float dimy, float dimz):m_Xbins(Xbins), m_Ybins(Ybins), m_Zbins(Zbins),m_xmin(xmin),m_xmax(xmax),m_ymin(ymin),m_ymax(ymax),m_zmin(zmin),m_zmax(zmax),m_dimx(dimx),m_dimy(dimy),m_dimz(dimz){}
  
  ~VoxelParam();

  float GetMinX(){return m_xmin;}
  float GetMinY(){return m_ymin;}
  float GetMinZ(){return m_zmin;}

  float GetMaxX(){return m_xmax;}
  float GetMaxY(){return m_ymax;}
  float GetMaxZ(){return m_zmax;}

  int GetZbins(){return m_Zbins;}
  float GetVoxelDimX(){return m_dimx;}
  float GetVoxelDimY(){return m_dimy;}
  float GetVoxelDimZ(){return m_dimz;}
  void CalculateFOV(GammaEvent *event);
  void FOVCalculator(std::vector<float> dir,std::vector<float> org, float theta, std::vector<float> &phiRanges, float rotAngY, float rotAngZ);
  float CalculatePhi(float x, float y);
 private:
  int m_Xbins, m_Ybins, m_Zbins;
  float m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax, m_dimx, m_dimy, m_dimz;

};
#endif
