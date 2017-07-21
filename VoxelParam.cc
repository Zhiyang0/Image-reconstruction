// Andrea Gutierrez
// University College London
// Medical Physics and Biomedical Engineering Department
// November 2016

#include "GammaEvent.hh"
#include "VoxelParam.hh"
#include "IntersectionParam.hh"
#include <iostream>
#include <list>
#define _USE_MATH_DEFINES
 
#include <cmath>
#define AXES 3
#define X_POS 0
#define Y_POS 1
#define Z_POS 2


#define XL 0
#define YT 1
#define XR 2
#define YB 3

#define EDGES 4

void VoxelParam::CalculateFOV(GammaEvent *event){
/* calculate min and max phi for max_z*/

  std::vector<float> phiRanges(2);
  //phiRanges[0] = phiMin
  //phiRanges[1] = phiMax

  FOVCalculator(event->GetDirectionVector(),event->GetOriginVector(),event->GetComptonAngle(), phiRanges, event->GetRotY(), event->GetRotZ());

  if(fabs(phiRanges[0]-phiRanges[1]) <0.03){// ranges too small
    phiRanges[0] = 0.0;
    phiRanges[1] = 2.0*M_PI;
  }
  
  event->SetPhiMin(phiRanges[0]);
  event->SetPhiMax(phiRanges[1]);
  
  //  std::cout << "small phi: "<< fabs(phiRanges[0]-phiRanges[1]) << " eventID: "<< event->GetEventID() << std::endl;

}

void VoxelParam::FOVCalculator(std::vector<float> dir,std::vector<float> org, float theta, std::vector<float> &phiRanges, float rotAngY, float rotAngZ){

//we need to solve quadratic equation equation discribing the intersection of the plane with the cone
//located at xlimit (S.J. Wilderman et al.)
//[nx(x-ax) + ny(y-ay) + nz(zlimit-az)]^2 = costheta^2[(x-ax)^2 + (y-ay)^2 + (zlimit -az)^2]
//where n is the unit vector along the cone axis
//ax, ay, az are the scatterer interaction position.
  float phiMin_zmin = 0.0;
  float phiMin_zmax = 0.0;
  
  float phiMax_zmin = 2*M_PI;
  float phiMax_zmax = 2*M_PI;
  float dirLength = sqrt(dir[0]*dir[0]  + dir[1]*dir[1] + dir[2]*dir[2]);
  
  std::vector<float> n(AXES);
  std::vector<int> nSolEdges_zmin(EDGES);
  std::vector<int> nSolEdges_zmax(EDGES);
  if (dirLength == 0) 
    dirLength = 1;
  
  for (int i = 0; i<AXES; i++){
    n[i]=dir[i]/dirLength;
  }
  
  
  // intersection class with quadratic solution for z_min and z_max:
  std::list<IntersectionParam*> intersect_zmin;
  std::list<IntersectionParam*> intersect_zmax;
  
  for (int i = 0; i<EDGES; i++){ // initialisation xmin, ymax, xmax, ymin
    IntersectionParam * interMin = new IntersectionParam(i);   
    intersect_zmin.push_back(interMin); 
    nSolEdges_zmin[i]=0;
    
    IntersectionParam * interMax = new IntersectionParam(i);   
    intersect_zmax.push_back(interMax);      
    nSolEdges_zmax[i]=0;
  }
  
   // solving four sides quadratic equations for zmin edge
  std::list<IntersectionParam*>::iterator itr_zmin;

  for(itr_zmin=intersect_zmin.begin(); itr_zmin!=intersect_zmin.end(); itr_zmin++)
    {
      (*itr_zmin)->CalculateQuadraticSolution(n,org,theta, m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, rotAngY, rotAngZ);
       nSolEdges_zmin[(*itr_zmin)->GetEdgeID()]=(*itr_zmin)->GetNSol();
    }
  

  // full circle in the plane, 0 <phi < 2*pi
  if (nSolEdges_zmin[XL]==0 &&  nSolEdges_zmin[XR]==0 && nSolEdges_zmin[YT]==0 && nSolEdges_zmin[YB]==0){
    phiMin_zmin = 0.0;
    phiMax_zmin = 2*M_PI;
  }
  
  //intersection in opposite sides XL and XR:
  else if (nSolEdges_zmin[XL]==1 &&  nSolEdges_zmin[XR]==1 && nSolEdges_zmin[YT]==0 && nSolEdges_zmin[YB]==0){
    //get intersection object
    itr_zmin = intersect_zmin.begin();
    std::advance(itr_zmin,XL);
    float temphi_1 = CalculatePhi((*itr_zmin)->GetSol1X(), (*itr_zmin)->GetSol1Y());
    itr_zmin = intersect_zmin.begin();
    std::advance(itr_zmin,XR);
    float temphi_2 = CalculatePhi((*itr_zmin)->GetSol1X(), (*itr_zmin)->GetSol1Y());
    //calculate phi from solution 
    if (temphi_1 < temphi_2){
      phiMin_zmin = temphi_1;
      phiMax_zmin = temphi_2;
    }
    else {
      phiMin_zmin = temphi_2;
      phiMax_zmin = temphi_1;
    }

    if(phiMax_zmin > (phiMin_zmin + M_PI)){
      temphi_1= phiMin_zmin;
      phiMin_zmin = phiMax_zmin -2*M_PI;
      phiMax_zmin = temphi_1;
    }

  }
  
  // intersection opppsite sides YT and YB
  else if (nSolEdges_zmin[XL]==0 &&  nSolEdges_zmin[XR]==0 && nSolEdges_zmin[YT]==1 && nSolEdges_zmin[YB]==1){
    //get intersection object
    itr_zmin = intersect_zmin.begin();
    std::advance(itr_zmin,YT);
    float temphi_1 = CalculatePhi((*itr_zmin)->GetSol1X(), (*itr_zmin)->GetSol1Y());
    itr_zmin = intersect_zmin.begin();
    std::advance(itr_zmin,YB);
    float temphi_2 = CalculatePhi((*itr_zmin)->GetSol1X(), (*itr_zmin)->GetSol1Y());
    //calculate phi from solution 
    if (temphi_1 < temphi_2){
      phiMin_zmin = temphi_1;
      phiMax_zmin = temphi_2;
    }
    else {
      phiMin_zmin = temphi_2;
      phiMax_zmin = temphi_1;
    }
  }


  // intersection adjacent sides, XL and YT
  else if (nSolEdges_zmin[XL]==1 &&  nSolEdges_zmin[XR]==0 && nSolEdges_zmin[YT]==1 && nSolEdges_zmin[YB]==0){
    //get intersection object
    itr_zmin = intersect_zmin.begin();
    std::advance(itr_zmin,XL);
    float temphi_1 = CalculatePhi((*itr_zmin)->GetSol1X(), (*itr_zmin)->GetSol1Y());
    itr_zmin = intersect_zmin.begin();
    std::advance(itr_zmin,YT);
    float temphi_2 = CalculatePhi((*itr_zmin)->GetSol1X(), (*itr_zmin)->GetSol1Y());
    //calculate phi from solution 
    if (temphi_1 < temphi_2){
      phiMin_zmin = temphi_1;
      phiMax_zmin = temphi_2;
    }
    else {
      phiMin_zmin = temphi_2;
      phiMax_zmin = temphi_1;
    }

    if(phiMax_zmin > (phiMin_zmin + M_PI)){
      temphi_1= phiMin_zmin;
      phiMin_zmin = phiMax_zmin -2*M_PI;
      phiMax_zmin = temphi_1;
    }
  }

 // intersection adjacent sides, XL and YB
  else if (nSolEdges_zmin[XL]==1 &&  nSolEdges_zmin[XR]==0 && nSolEdges_zmin[YT]==0 && nSolEdges_zmin[YB]==1){
    //get intersection object
    itr_zmin = intersect_zmin.begin();
    std::advance(itr_zmin,XL);
    float temphi_1 = CalculatePhi((*itr_zmin)->GetSol1X(), (*itr_zmin)->GetSol1Y());
    itr_zmin = intersect_zmin.begin();
    std::advance(itr_zmin,YB);
    float temphi_2 = CalculatePhi((*itr_zmin)->GetSol1X(), (*itr_zmin)->GetSol1Y());
    //calculate phi from solution 
    if (temphi_1 < temphi_2){
      phiMin_zmin = temphi_1;
      phiMax_zmin = temphi_2;
    }
    else {
      phiMin_zmin = temphi_2;
      phiMax_zmin = temphi_1;
    }
    if(phiMax_zmin > (phiMin_zmin + M_PI)){
      temphi_1= phiMin_zmin;
      phiMin_zmin = phiMax_zmin -2*M_PI;
      phiMax_zmin = temphi_1;
    }
  }


 // intersection adjacent sides, XR and YT
  else if (nSolEdges_zmin[XL]==0 &&  nSolEdges_zmin[XR]==1 && nSolEdges_zmin[YT]==1 && nSolEdges_zmin[YB]==0){
    //get intersection object
    itr_zmin = intersect_zmin.begin();
    std::advance(itr_zmin,XR);
    float temphi_1 = CalculatePhi((*itr_zmin)->GetSol1X(), (*itr_zmin)->GetSol1Y());
    itr_zmin = intersect_zmin.begin();
    std::advance(itr_zmin,YT);
    float temphi_2 = CalculatePhi((*itr_zmin)->GetSol1X(), (*itr_zmin)->GetSol1Y());
    //calculate phi from solution 
    if (temphi_1 < temphi_2){
      phiMin_zmin = temphi_1;
      phiMax_zmin = temphi_2;
    }
    else {
      phiMin_zmin = temphi_2;
      phiMax_zmin = temphi_1;
    }

    if(phiMax_zmin > (phiMin_zmin + M_PI)){
      temphi_1= phiMin_zmin;
      phiMin_zmin = phiMax_zmin -2*M_PI;
      phiMax_zmin = temphi_1;
    }
  }

 // intersection adjacent sides, XR and YB
  else if (nSolEdges_zmin[XL]==0 &&  nSolEdges_zmin[XR]==1 && nSolEdges_zmin[YT]==0 && nSolEdges_zmin[YB]==1){
    //get intersection object
    itr_zmin = intersect_zmin.begin();
    std::advance(itr_zmin,XR);
    float temphi_1 = CalculatePhi((*itr_zmin)->GetSol1X(), (*itr_zmin)->GetSol1Y());
    itr_zmin = intersect_zmin.begin();
    std::advance(itr_zmin,YB);
    float temphi_2 = CalculatePhi((*itr_zmin)->GetSol1X(), (*itr_zmin)->GetSol1Y());
    //calculate phi from solution 
    if (temphi_1 < temphi_2){
      phiMin_zmin = temphi_1;
      phiMax_zmin = temphi_2;
    }
    else {
      phiMin_zmin = temphi_2;
      phiMax_zmin = temphi_1;
    }

    if(phiMax_zmin > (phiMin_zmin + M_PI)){
      temphi_1= phiMin_zmin;
      phiMin_zmin = phiMax_zmin -2*M_PI;
      phiMax_zmin = temphi_1;
    }
  }




  ////////////////////////////////////////////////////
  // solving four sides quadratic equations for zmax edge
  std::list<IntersectionParam*>::iterator itr_zmax;
   for(itr_zmax=intersect_zmax.begin(); itr_zmax!=intersect_zmax.end(); itr_zmax++)
    {
      (*itr_zmax)->CalculateQuadraticSolution(n,org,theta, m_xmin, m_xmax, m_ymin, m_ymax, m_zmax, rotAngY, rotAngZ);
       nSolEdges_zmax[(*itr_zmax)->GetEdgeID()]=(*itr_zmax)->GetNSol();
    }
  

  // full circle in the plane, 0 <phi < 2*pi
  if (nSolEdges_zmax[XL]==0 &&  nSolEdges_zmax[XR]==0 && nSolEdges_zmax[YT]==0 && nSolEdges_zmax[YB]==0){
    phiMin_zmax = 0.0;
    phiMax_zmax = 2*M_PI;
  }



  //intersection in opposite sides XL and XR:
  else if (nSolEdges_zmax[XL]==1 &&  nSolEdges_zmax[XR]==1 && nSolEdges_zmax[YT]==0 && nSolEdges_zmax[YB]==0){
    //get intersection object
    itr_zmax = intersect_zmax.begin();
    std::advance(itr_zmax,XL);
    float temphi_1 = CalculatePhi((*itr_zmax)->GetSol1X(), (*itr_zmax)->GetSol1Y());
    itr_zmax = intersect_zmax.begin();
    std::advance(itr_zmax,XR);
    float temphi_2 = CalculatePhi((*itr_zmax)->GetSol1X(), (*itr_zmax)->GetSol1Y());
    //calculate phi from solution 
    if (temphi_1 < temphi_2){
      phiMin_zmax = temphi_1;
      phiMax_zmax = temphi_2;
    }
    else {
      phiMin_zmax = temphi_2;
      phiMax_zmax = temphi_1;
    }

    if(phiMax_zmax > (phiMin_zmax + M_PI)){
      temphi_1= phiMin_zmax;
      phiMin_zmax = phiMax_zmax -2*M_PI;
      phiMax_zmax = temphi_1;
    }
  }
  
  // intersection opppsite sides YT and YB
  else if (nSolEdges_zmax[XL]==0 &&  nSolEdges_zmax[XR]==0 && nSolEdges_zmax[YT]==1 && nSolEdges_zmax[YB]==1){
    //get intersection object
    itr_zmax = intersect_zmax.begin();
    std::advance(itr_zmax,YT);
    float temphi_1 = CalculatePhi((*itr_zmax)->GetSol1X(), (*itr_zmax)->GetSol1Y());
    itr_zmax = intersect_zmax.begin();
    std::advance(itr_zmax,YB);
    float temphi_2 = CalculatePhi((*itr_zmax)->GetSol1X(), (*itr_zmax)->GetSol1Y());
    //calculate phi from solution 
    if (temphi_1 < temphi_2){
      phiMin_zmax = temphi_1;
      phiMax_zmax = temphi_2;
    }
    else {
      phiMin_zmax = temphi_2;
      phiMax_zmax = temphi_1;
    }
    if(phiMax_zmax > (phiMin_zmax + M_PI)){
      temphi_1= phiMin_zmax;
      phiMin_zmax = phiMax_zmax -2*M_PI;
      phiMax_zmax = temphi_1;
    }
  }


  // intersection adjacent sides, XL and YT
  else if (nSolEdges_zmax[XL]==1 &&  nSolEdges_zmax[XR]==0 && nSolEdges_zmax[YT]==1 && nSolEdges_zmax[YB]==0){
    //get intersection object
    itr_zmax = intersect_zmax.begin();
    std::advance(itr_zmax,XL);
    float temphi_1 = CalculatePhi((*itr_zmax)->GetSol1X(), (*itr_zmax)->GetSol1Y());
    itr_zmax = intersect_zmax.begin();
    std::advance(itr_zmax,YT);
    float temphi_2 = CalculatePhi((*itr_zmax)->GetSol1X(), (*itr_zmax)->GetSol1Y());
    //calculate phi from solution 
    if (temphi_1 < temphi_2){
      phiMin_zmax = temphi_1;
      phiMax_zmax = temphi_2;
    }
    else {
      phiMin_zmax = temphi_2;
      phiMax_zmax = temphi_1;
    }

    if(phiMax_zmax > (phiMin_zmax + M_PI)){
      temphi_1= phiMin_zmax;
      phiMin_zmax = phiMax_zmax -2*M_PI;
      phiMax_zmax = temphi_1;
    }
  }

 // intersection adjacent sides, XL and YB
  else if (nSolEdges_zmax[XL]==1 &&  nSolEdges_zmax[XR]==0 && nSolEdges_zmax[YT]==0 && nSolEdges_zmax[YB]==1){
    //get intersection object
    itr_zmax = intersect_zmax.begin();
    std::advance(itr_zmax,XL);
    float temphi_1 = CalculatePhi((*itr_zmax)->GetSol1X(), (*itr_zmax)->GetSol1Y());
    itr_zmax = intersect_zmax.begin();
    std::advance(itr_zmax,YB);
    float temphi_2 = CalculatePhi((*itr_zmax)->GetSol1X(), (*itr_zmax)->GetSol1Y());
    //calculate phi from solution 
    if (temphi_1 < temphi_2){
      phiMin_zmax = temphi_1;
      phiMax_zmax = temphi_2;
    }
    else {
      phiMin_zmax = temphi_2;
      phiMax_zmax = temphi_1;
    }

    if(phiMax_zmax > (phiMin_zmax + M_PI)){
      temphi_1= phiMin_zmax;
      phiMin_zmax = phiMax_zmax -2*M_PI;
      phiMax_zmax = temphi_1;
    }
  }


 // intersection adjacent sides, XR and YT
  else if (nSolEdges_zmax[XL]==0 &&  nSolEdges_zmax[XR]==1 && nSolEdges_zmax[YT]==1 && nSolEdges_zmax[YB]==0){
    //get intersection object
    itr_zmax = intersect_zmax.begin();
    std::advance(itr_zmax,XR);
    float temphi_1 = CalculatePhi((*itr_zmax)->GetSol1X(), (*itr_zmax)->GetSol1Y());
    itr_zmax = intersect_zmax.begin();
    std::advance(itr_zmax,YT);
    float temphi_2 = CalculatePhi((*itr_zmax)->GetSol1X(), (*itr_zmax)->GetSol1Y());
    //calculate phi from solution 
    if (temphi_1 < temphi_2){
      phiMin_zmax = temphi_1;
      phiMax_zmax = temphi_2;
    }
    else {
      phiMin_zmax = temphi_2;
      phiMax_zmax = temphi_1;
    }

    if(phiMax_zmax > (phiMin_zmax + M_PI)){
      temphi_1= phiMin_zmax;
      phiMin_zmax = phiMax_zmax -2*M_PI;
      phiMax_zmax = temphi_1;
    }
  }

 // intersection adjacent sides, XR and YB
  else if (nSolEdges_zmax[XL]==0 &&  nSolEdges_zmax[XR]==1 && nSolEdges_zmax[YT]==0 && nSolEdges_zmax[YB]==1){
    //get intersection object
    itr_zmax = intersect_zmax.begin();
    std::advance(itr_zmax,XR);
    float temphi_1 = CalculatePhi((*itr_zmax)->GetSol1X(), (*itr_zmax)->GetSol1Y());
    itr_zmax = intersect_zmax.begin();
    std::advance(itr_zmax,YB);
    float temphi_2 = CalculatePhi((*itr_zmax)->GetSol1X(), (*itr_zmax)->GetSol1Y());
    //calculate phi from solution 
    if (temphi_1 < temphi_2){
      phiMin_zmax = temphi_1;
      phiMax_zmax = temphi_2;
    }
    else {
      phiMin_zmax = temphi_2;
      phiMax_zmax = temphi_1;
    }

    if(phiMax_zmax > (phiMin_zmax + M_PI)){
      temphi_1= phiMin_zmax;
      phiMin_zmax = phiMax_zmax -2*M_PI;
      phiMax_zmax = temphi_1;
    }
  }

  if(phiMin_zmin < phiMin_zmax){
    phiRanges[0]=phiMin_zmin;
  }
  else phiRanges[0]=phiMin_zmax;

  if(phiMax_zmin > phiMax_zmax){
    phiRanges[1]=phiMax_zmin;
  }
  else phiRanges[1]=phiMax_zmax;


}

float VoxelParam::CalculatePhi(float x, float y){
  float phi;
  if(x == 0)
    {
      //y is positive --> phi = pi/2
      if(y > 0)
	phi = M_PI/2.0;
	
      else// y is negative =--> 3*phi/2
	phi = 3*M_PI/2.0;
    }
  else if(y == 0){
    if (x > 0)
      phi = 0.0;
    else 
      phi = M_PI;
  }
  else {
    float tanphi = fabs(y/x);// tan(phi) = opp/adj
    if( x > 0 && y > 0){
      phi = atan(tanphi);
    }
    if(x < 0 && y > 0){
      phi = M_PI - atan(tanphi);
    }
    if( x < 0 && y < 0){
      phi = M_PI +  atan(tanphi);
    }
    if( x > 0 && y < 0){
      phi = 2*M_PI - atan(tanphi);
    }
  }
  return phi;
}
