// Andrea Gutierrez
// University College London
// Medical Physics and Biomedical Engineering Department
// November 2016

#include "IntersectionParam.hh"
#include <math.h>

#define XL 0
#define YT 1
#define XR 2
#define YB 3

#define X_POS 0
#define Y_POS 1
#define Z_POS 2


void IntersectionParam::CalculateQuadraticSolution(std::vector<float> &n, std::vector<float> &org, float theta, float xmin, float xmax, float ymin, float ymax, float z, float rotAngY, float rotAngZ){

  float leftTerm, rightTerm, aCoeff, bCoeff, cCoeff;
  
  if(edgeID == XL){ //solve quadratic equation for X left (xmin)
    
    leftTerm = CalculateLeftTerm(n[X_POS],xmin, org[X_POS], n[Z_POS],z,org[Z_POS]);
    rightTerm = CalculateRightTerm(n[X_POS],xmin, org[X_POS], n[Z_POS],z,org[Z_POS]); 
    aCoeff = CalculateACoeff(n[Y_POS], theta);
    bCoeff = CalculateBCoeff(n[Y_POS], leftTerm);
    cCoeff = CalculateCCoeff(leftTerm, rightTerm, theta);
    SolveQuadraticEq(aCoeff, bCoeff, cCoeff);
    InLimitsSolutions(ymin, ymax, org[Y_POS]);
    SetSolutionVector(xmin,z,true, rotAngY, rotAngZ, org, theta);//limitX = true
  }
  else if(edgeID == YT){
    leftTerm = CalculateLeftTerm(n[Y_POS],ymax, org[Y_POS], n[Z_POS],z,org[Z_POS]);
    rightTerm = CalculateRightTerm(n[Y_POS],ymax, org[Y_POS], n[Z_POS],z,org[Z_POS]); 
    aCoeff = CalculateACoeff(n[X_POS], theta);
    bCoeff = CalculateBCoeff(n[X_POS], leftTerm);
    cCoeff = CalculateCCoeff(leftTerm, rightTerm, theta);
    SolveQuadraticEq(aCoeff, bCoeff, cCoeff);
    InLimitsSolutions(xmin, xmax, org[X_POS]);
    SetSolutionVector(ymax,z,false, rotAngY, rotAngZ, org, theta );//limitX = false
  }
  else if(edgeID == XR){ //solve quadratic equation for X left (xmax)
    
    leftTerm = CalculateLeftTerm(n[X_POS],xmax, org[X_POS], n[Z_POS],z,org[Z_POS]);
    rightTerm = CalculateRightTerm(n[X_POS],xmax, org[X_POS], n[Z_POS],z,org[Z_POS]); 
    aCoeff = CalculateACoeff(n[Y_POS], theta);
    bCoeff = CalculateBCoeff(n[Y_POS], leftTerm);
    cCoeff = CalculateCCoeff(leftTerm, rightTerm, theta);
    SolveQuadraticEq(aCoeff, bCoeff, cCoeff);
    InLimitsSolutions(ymin, ymax, org[Y_POS]);
    SetSolutionVector(xmax,z,true, rotAngY, rotAngZ, org, theta);//limitX = true
  }
  else if(edgeID == YB){
    leftTerm = CalculateLeftTerm(n[Y_POS],ymin, org[Y_POS], n[Z_POS],z,org[Z_POS]);
    rightTerm = CalculateRightTerm(n[Y_POS],ymin, org[Y_POS], n[Z_POS],z,org[Z_POS]); 
    aCoeff = CalculateACoeff(n[X_POS], theta);
    bCoeff = CalculateBCoeff(n[X_POS], leftTerm);
    cCoeff = CalculateCCoeff(leftTerm, rightTerm, theta);
    SolveQuadraticEq(aCoeff, bCoeff, cCoeff);
    InLimitsSolutions(xmin, xmax, org[X_POS]);
    SetSolutionVector(ymin,z,false, rotAngY, rotAngZ, org, theta);//limitX = false
  }

  
  //works up to here;



}

  float IntersectionParam::CalculateLeftTerm(float n_j,float limit_j, float org_j, float n_z, float limit_z, float org_z){
  
  float leftTerm;
  leftTerm = n_j*(limit_j-org_j) + n_z*(limit_z - org_z);
  return leftTerm; 

}

  float IntersectionParam::CalculateRightTerm(float n_j,float limit_j, float org_j, float n_z, float limit_z, float org_z){
  
  float rightTerm;
  rightTerm = (limit_j-org_j)*(limit_j-org_j) + (limit_z - org_z)*(limit_z - org_z);
  return rightTerm; 

  }

float IntersectionParam::CalculateACoeff(float n_k, float theta){

  return ((n_k*n_k) - (cos(theta)*cos(theta)));
}


float IntersectionParam::CalculateBCoeff(float n_k, float l){
  return (2*l*n_k);
}

float IntersectionParam::CalculateCCoeff(float l,float r, float theta){
  return (l*l - cos(theta)*cos(theta)*r);

}

void IntersectionParam::SolveQuadraticEq(float a, float b, float c){

  float disc = (b*b) - (4*a*c);
  
  if (disc>= 0){
    q1 = (-1*b +sqrt(disc))/(2*a); //first solution or unque if disc == 0
    nQSol++;
  }
  if(disc > 0){
    q2 = (-1*b -sqrt(disc))/(2*a); //second solution
    nQSol++;
  }

}

void IntersectionParam::InLimitsSolutions(float min, float max, float org){

   if (nQSol > 0 && nSol == 0){
      if (q1 > min && q1 < max){
	nSol++;
	s1 = q1 + org;
      }
    }


    if (nQSol == 2 && nSol == 1){
      if (q2 > min && q2 < max){
	nSol++;
	s2 = q2 + org;
      } 
    }
    
    
    
    if (nQSol == 2 && nSol == 0){
      if (q2 > min && q2 < max){
	nSol++;
	s1 = q2 + org;
      } 
    }
}


void IntersectionParam::SetSolutionVector(float limit, float z, bool limitX, float rotAngY, float rotAngZ, std::vector<float> &org, float theta){
 
  std::vector<float> tmp1(AXES);
  std::vector<float> tmp2(AXES);

  std::vector<float> unrotZ(AXES);



  if(limitX){
    if(nSol > 0){
      tmp1[0] = limit - org[X_POS];
      tmp1[1] = s1 - org[Y_POS];
      tmp1[2] = z - org[Z_POS];
    }
    if(nSol ==2){
      tmp2[0] = limit - org[X_POS];
      tmp2[1] = s2 - org[Y_POS];
      tmp2[2] = z - org[Z_POS];
    }
  }
  else {
    if(nSol > 0){
      tmp1[0] = s1 - org[X_POS];
      tmp1[1] = limit - org[Y_POS];
      tmp1[2] = z - org[Z_POS];
    }
    if(nSol ==2){
      tmp2[0] = s2 - org[X_POS];
      tmp2[1] = limit - org[Y_POS];
      tmp2[2] = z - org[Z_POS];
    }
  }

  
  std::vector<std::vector<float> > rotZ(AXES, std::vector<float>(AXES, 0));
  std::vector<std::vector<float> > rotY(AXES, std::vector<float>(AXES, 0));
  // inversed rotation matrix on Y
  rotY[0][0] = cos(-rotAngY);
  rotY[1][0] = 0;
  rotY[2][0] = sin(-rotAngY);
  
  rotY[0][1] = 0;
  rotY[1][1] = 1;
  rotY[2][1] = 0;
  
  rotY[0][2] = -sin(-rotAngY);
  rotY[1][2] = 0.0;
  rotY[2][2] = cos(-rotAngY);
  
  //inversed rotation matrix on Z axis 
 rotZ[0][0] = cos(-rotAngZ);
 rotZ[1][0] =-sin(-rotAngZ);
 rotZ[2][0] = 0;
 
 rotZ[0][1] = sin(-rotAngZ);
 rotZ[1][1] = cos(-rotAngZ);
 rotZ[2][1] = 0;

 rotZ[0][2] = 0;
 rotZ[1][2] = 0.0;
 rotZ[2][2] = 1;

 if(nSol>0){ 
  
   unrotZ[0]=rotZ[0][0]*tmp1[0] + rotZ[1][0]*tmp1[1]+ + rotZ[2][0]*tmp1[2];
   unrotZ[1]=rotZ[0][1]*tmp1[0] + rotZ[1][1]*tmp1[1]+ + rotZ[2][1]*tmp1[2];
   unrotZ[2]=rotZ[0][2]*tmp1[0] + rotZ[1][2]*tmp1[1]+ + rotZ[2][2]*tmp1[2];    
   
   
   tmp1[0]=rotY[0][0]*unrotZ[0] + rotY[1][0]*unrotZ[1]+ + rotY[2][0]*unrotZ[2];
   tmp1[1]=rotY[0][1]*unrotZ[0] + rotY[1][1]*unrotZ[1]+ + rotY[2][1]*unrotZ[2];
   tmp1[2]=rotY[0][2]*unrotZ[0] + rotY[1][2]*unrotZ[1]+ + rotY[2][2]*unrotZ[2];   

   if(tmp1[Z_POS] != 0){
     sol1[0] = tmp1[0]*cos(theta)/tmp1[Z_POS];
     sol1[1] = tmp1[1]*cos(theta)/tmp1[Z_POS];
     sol1[2] = tmp1[2]*cos(theta)/tmp1[Z_POS];
   }
  
 }

 if(nSol ==2){

   unrotZ[0]=rotZ[0][0]*tmp2[0] + rotZ[1][0]*tmp2[1]+ + rotZ[2][0]*tmp2[2];
   unrotZ[1]=rotZ[0][1]*tmp2[0] + rotZ[1][1]*tmp2[1]+ + rotZ[2][1]*tmp2[2];
   unrotZ[2]=rotZ[0][2]*tmp2[0] + rotZ[1][2]*tmp2[1]+ + rotZ[2][2]*tmp2[2];    
   
   
   tmp2[0]=rotY[0][0]*unrotZ[0] + rotY[1][0]*unrotZ[1]+ + rotY[2][0]*unrotZ[2];
   tmp2[1]=rotY[0][1]*unrotZ[0] + rotY[1][1]*unrotZ[1]+ + rotY[2][1]*unrotZ[2];
   tmp2[2]=rotY[0][2]*unrotZ[0] + rotY[1][2]*unrotZ[1]+ + rotY[2][2]*unrotZ[2];  

    if(tmp2[Z_POS] != 0){
     sol2[0] = tmp2[0]*cos(theta)/tmp2[Z_POS];
     sol2[1] = tmp2[1]*cos(theta)/tmp2[Z_POS];
     sol2[2] = tmp2[2]*cos(theta)/tmp2[Z_POS];
   }

  
 }


 

}

