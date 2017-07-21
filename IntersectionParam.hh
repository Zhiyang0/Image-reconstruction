// Andrea Gutierrez
// University College London
// Medical Physics and Biomedical Engineering Department
// November 2016

#ifndef INTERSECTIONPARAM_H
#define INTERSECTIONPARAM_H
#include <vector>
#include <iostream>
#define AXES 3

class IntersectionParam
{
public:
  IntersectionParam(int edge):edgeID(edge),sol1(AXES),sol2(AXES){
    nSol = 0;
    
    s1 = -1;
    s2 = -1;

    sol1[0]= -1;
    sol1[1]= -1;
    sol1[2]= -1;
    
    sol2[0]= -1;
    sol2[1]= -1;
    sol2[2]= -1;
    
    phi_min = -1;
    phi_max = -1;
    
    q1 = -1;
    q2 = -1;
  
    nQSol = 0;
  }
  ~IntersectionParam(){};

  int GetNSol(){return nSol;}
  int GetNQSol(){return nQSol;}

  void SetNSol(int sol){nSol = sol;}
  int GetEdgeID(){return edgeID;}

  float GetSol1X(){return sol1[0];}
  float GetSol1Y(){return sol1[1];}

  float GetSol2X(){return sol2[0];}
  float GetSol2Y(){return sol2[1];}


  void CalculateQuadraticSolution(std::vector<float> &n, std::vector<float> &org, float theta, float xmin, float xmax, float ymin, float ymax, float z, float rotAngY, float rotAngZ);
  float CalculateLeftTerm(float n_j,float limit_j, float org_j, float n_z, float limit_z, float org_z);
  float CalculateRightTerm(float n_j,float limit_j, float org_j, float n_z, float limit_z, float org_z);
  float CalculateACoeff(float n_k, float theta);
  float CalculateBCoeff(float n_k, float l);
  float CalculateCCoeff(float l, float r, float theta);
  void SolveQuadraticEq(float a, float b, float c);
  void InLimitsSolutions(float min, float max, float org);
  void SetSolutionVector(float limit, float z, bool limitX, float rotAngY, float rotAngZ, std::vector<float> &org, float theta);
  
private:
  //solutions to be used to find phi ranges
  int nSol;
  int edgeID;
  std::vector<float> sol1;
  std::vector<float> sol2;
  
  float s1, s2;
  float phi_min;
  float phi_max; 

  //results from quadratic equation
  float q1;
  float q2;
  int nQSol;

};
#endif
