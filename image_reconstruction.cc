// Andrea Gutierrez
// University College London
// Medical Physics and Biomedical Engineering Department
// November 2016

//SOE algorithm: Andreyev A, Sitek A, Celler A, Fast image reconstruction for Compton camera using stochastic origin ensemble approach. Med Phys. 2011 Jan;38(1):429-38.
//Code implementation: Andrea Gutierrez, UCL

// test by xiao in THU

#include <iostream>
#include <fstream>
#include <list>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "GammaEvent.hh"

#include <TRandom1.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TNtuple.h>

#define X_POS 0
#define Y_POS 1
#define Z_POS 2

#define ELECTRON_MASS 510.999 // electron mass in keV
#define AXES 3
#define SOE_ITR 1
using namespace std;

class GammaEvent;
class Mat
{
public:
  Mat(int x1,int y1,int z1);
  Mat();
  ~Mat() {cout<<"DESTRUCT\n";};
  void Matrix();
  int& at(int x1,int y1,int z1);
  int *d;
private:
  int x,y,z;
  int size;
};
Mat::Mat()
{}
Mat::Mat(int x1,int y1,int z1):x(x1),y(y1),z(z1)
{
  size=x*y*z;
  d = new int[size]();
}
int &Mat::at(int x_,int y_,int z_)
{
  int index=z_+z*(y_+y*x_);
  assert(index<size);
  return d[index];
}

bool ComptonEdgeTest (float scatt_e1,float abs_e2)
{
  bool pass;
  //calculate compton edge test:
  float E_i = (scatt_e1+abs_e2);
  float E_max = 2*E_i*E_i/(ELECTRON_MASS + 2*E_i);
  if(scatt_e1 > E_max)   
    pass = false;
  else pass = true;
  return pass;
}  


int main()
{
  clock_t t_start, t_init, t_itr;
  t_start  = clock();
  ifstream infile;

  infile.open("simData.txt",ios::in); //this one to be used
 
  int Xbins = 400;//1000;
  int Ybins = 400;//1000;
  int Zbins = 400;
  float X_limit=Xbins;
  float Y_limit=Ybins*cos(0.167*3.1415);
  float Z_limit=Zbins*sin(0.167*3.1415);
  

  float xmin = -200;// mm
  float xmax = 200; // mm
  float ymin = -200;// mm
  float ymax = 200; // mm
  float zmin = -250; // mm
  float zmax = 150;// mm
  
  float voxel_dimX =(xmax-xmin)/float(Xbins); // 1 mm.
  float voxel_dimY =(ymax-ymin)/float(Ybins); // 1 mm.
  float voxel_dimZ =(zmax-zmin)/float(Zbins); // 1 mm.
 
  VoxelParam *imageParam = new VoxelParam(Xbins,Ybins,Zbins,xmin,xmax,ymin,ymax,zmin,zmax,voxel_dimX,voxel_dimY,voxel_dimZ);

  std::cout << "Get min image: " <<  imageParam->GetMinX() << std::endl;

  std::cout << "Get voxel dim image: " <<  imageParam->GetVoxelDimX() << std::endl;\

  Mat pDensity(Xbins,Ybins,Zbins);
  Mat finalImage(Xbins,Ybins,Zbins);
 // float ***pDensity;
 // pDensity = new float**[Zbins];

 // float ***finalImage;
 // finalImage = new float**[Zbins];

 // for (int zitr =0;zitr < Zbins;zitr++){
 //   pDensity[zitr] = new float*[Ybins];
 //   finalImage[zitr] = new float*[Ybins];
 //   for(int yitr = 0; yitr < Ybins; yitr++){
 //     pDensity[zitr][yitr] = new float[Xbins];
 //     finalImage[zitr][yitr] = new float[Xbins];
 //   }
 // }

 // //initialise density matrix to 0;
 // for (int xitr =0;xitr < Xbins;xitr++)
 //   for (int yitr =0;yitr < Ybins;yitr++)
 //     for (int zitr =0;zitr < Zbins;zitr++){
 //       pDensity[zitr][yitr][xitr]=0;
 //       finalImage[zitr][yitr][xitr]=0;
 //     }
  
  float x1,y1,z1,e1,x2,y2,z2,e2;
  float eventItr = 0;
  int failed = 0;
  int outMatrix = 0;
  std::list<GammaEvent*> gammaSet;
  bool passCEtest = false; // Compton Edge Test
    
  //while (!infile.eof() && eventItr<9){
  while (!infile.eof()){

    infile >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> e1 >> e2;


    //-- IMPORTANT ----- we have considered that the data is given with right interaction sequence
    //if the event pass the compton edge test, then the event is stored.
    //E_edge = 2*E_gamma/(electron_mass +2*E_gammma).
    //This is the maximum energy a photon can deposit in a compton scattering event.
    
    passCEtest = ComptonEdgeTest(e1,e2);    
    if ((passCEtest)){  //116 126 and 160, 170
      //if ((passCEtest) && ((e1+e2)>160 && (e1+e2)<170)){  //116 126 and 160, 170
    //if ((passCEtest) && (((e1+e2)>160 && (e1+e2)<170) ||  ((e1+e2)>116 && (e1+e2)<126))){  //116 126 and 160, 170  //116 126 and 160, 170

       //    for(int k = 0;k<1;k++){//500 //200
      //using 2 sources
      GammaEvent * event= new GammaEvent(x1,y1,z1,e1,x2,y2,z2,e2);
      gammaSet.push_back(event);
      eventItr++;
      // }
      
    } 
  }
  std::cout << "Total number of events: " << eventItr << std::endl;
  infile.close();
  std::cout << "step 1, reading events:  " << float((clock() - t_start))/CLOCKS_PER_SEC << " sec" << std::endl;
  t_init=clock();
 

 std::vector<float> Z_sum(Zbins);
  for (int zitr =0;zitr < Zbins;zitr++){
    Z_sum[zitr] = 0;
  }


  std::vector<int> newVoxelPos(AXES);
  std::vector<int> oldVoxelPos(AXES);
  std::vector<float> new_pos(AXES);
  std::vector<float> old_pos(AXES);

  std::list<GammaEvent*>::iterator itr;
  /*initial iteration over all the gamma events to calculate:
    - compton angle
    - conical surface rotation, origin, etc
    - new random position
   */
  float event = 1;

  for(itr=gammaSet.begin(); itr!=gammaSet.end(); itr++)
    {
      //set event ID
      (*itr)->SetEventID(event);
      // calculate compton angle.
      (*itr)->CalculateComptonAngle();
      //calculate conical surface parameters only once e.g. rotation matrix
      (*itr)->CalculateConicalSurface();
      imageParam->CalculateFOV((*itr));

    
      for (int k = 0;k<1; k++){//500
      //uses random numbersfrom stdlib     
      (*itr)->CalculateNewRandomPosition(new_pos,imageParam);

      // set positions
      (*itr)->SetOldPos(new_pos);      
      (*itr)->SetNewPos(new_pos);
      
      //find voxel position -- use fonction instead?
      newVoxelPos[X_POS] = int((new_pos[X_POS] - xmin)/voxel_dimX);
      newVoxelPos[Y_POS] = int((new_pos[Y_POS] - ymin)/voxel_dimY);
      newVoxelPos[Z_POS] = int((new_pos[Z_POS] - zmin)/voxel_dimZ);

    
      if ((newVoxelPos[X_POS] <= 0)||(newVoxelPos[Y_POS] <= 0)||(newVoxelPos[Z_POS] < 0)||(newVoxelPos[X_POS] >= Xbins) || (newVoxelPos[Y_POS]>=Ybins) || (newVoxelPos[Z_POS]>=Zbins)){
	outMatrix++;
	//std::cout << "oustide range " << std::endl;
	}

	else {
	  pDensity.at(newVoxelPos[X_POS],newVoxelPos[Y_POS],newVoxelPos[Z_POS])++;

	  }
      //std::cout << "\n" << std::endl;

    }

      event++;
    }//end of initilisation over Gamma Data Set.

  std::cout << "Read " << event << " Compton events" << std::endl;
  std::cout << failed << " events failed Compton edge test" << std::endl;
  std::cout << "step 2, cone calculation:  " << float((clock() - t_init))/CLOCKS_PER_SEC << " sec" <<std::endl;



  /****** checking Z *****/
  /* for (int zitr =0;zitr < Zbins;zitr++){
    for (int yitr =0;yitr < Ybins;yitr++){
      for (int xitr =0;xitr < Xbins;xitr++){
	Z_sum[zitr]+=pDensity[zitr][yitr][zitr];
      }
    }
  }
  */
  std::cout << "Reading objects from list " << std::endl;
  std::cout << outMatrix <<"\t events out of the image space"  << std::endl;
  
  float newProb,oldProb, accept;
  bool oldOutside, newOutside, isAccepted;
  float randomNum;
  
 // main iteration of the SOE algori0tm
  for(int algo_itr = 0; algo_itr <SOE_ITR; algo_itr++){
    //  std::cout << "Starting iteration: " << algo_itr << std::endl;
    for(itr=gammaSet.begin(); itr!=gammaSet.end(); itr++)
      {	
	//get old positions/voxels:
	old_pos[X_POS]=(*itr)->GetOldX();
	old_pos[Y_POS]=(*itr)->GetOldY();
	old_pos[Z_POS]=(*itr)->GetOldZ();

	oldVoxelPos[X_POS] = int((old_pos[X_POS] - xmin)/voxel_dimX);
	oldVoxelPos[Y_POS] = int((old_pos[Y_POS] - ymin)/voxel_dimY);
	oldVoxelPos[Z_POS] = int((old_pos[Z_POS] - zmin)/voxel_dimZ);
	//for(int k=0;k<1000;k++){
	
	//calculate and set new random position
	(*itr)->CalculateNewRandomPosition(new_pos, imageParam);
	(*itr)->SetNewPos(new_pos);

	//get voxel position of new position
	newVoxelPos[X_POS] = int((new_pos[X_POS] - xmin)/voxel_dimX);
	newVoxelPos[Y_POS] = int((new_pos[Y_POS] - ymin)/voxel_dimY);
	newVoxelPos[Z_POS] = int((new_pos[Z_POS] - zmin)/voxel_dimZ);
	isAccepted = false;
	if ((newVoxelPos[X_POS] <= 0)||(newVoxelPos[Y_POS] <= 0)||(newVoxelPos[Z_POS] < 0)||(newVoxelPos[X_POS] >= Xbins) || (newVoxelPos[Y_POS]>=Ybins) || (newVoxelPos[Z_POS]>=Zbins)){
	  newProb=0;
	  newOutside = true;
	}
	else {
	  newProb =  pDensity.at(newVoxelPos[X_POS],newVoxelPos[Y_POS],newVoxelPos[Z_POS])+1;
	  newOutside = false;
	  
	}

	if ((oldVoxelPos[X_POS] <= 0)||(oldVoxelPos[Y_POS] <= 0)||(oldVoxelPos[Z_POS] < 0)||(oldVoxelPos[X_POS] >= Xbins) || (oldVoxelPos[Y_POS]>=Ybins) || (oldVoxelPos[Z_POS]>=Zbins)){
	  oldProb=0;
	  oldOutside = true;
	}
	else {
	  oldProb =  pDensity.at(oldVoxelPos[X_POS],oldVoxelPos[Y_POS],oldVoxelPos[Z_POS]);
	  oldOutside = false;
	}

	if (oldProb>0){
	  accept =newProb/oldProb;
	}
	else {accept = 1;}

	//accept if higher than or equal 1
	if (accept >=1 && newOutside == false){
	  isAccepted = true;
	  //increase density at the new position
	  pDensity.at(newVoxelPos[X_POS],newVoxelPos[Y_POS],newVoxelPos[Z_POS])++;
	  Z_sum[newVoxelPos[Z_POS]]++;
	  //new position becomes old position 
	  (*itr)->SetOldPos(new_pos);

	  if(oldOutside==false){
	  //decrease density at old position
	    pDensity.at(oldVoxelPos[X_POS],oldVoxelPos[Y_POS],oldVoxelPos[Z_POS])--;
	    Z_sum[newVoxelPos[Z_POS]]--;
	  }
	}
	if(accept <1 && newOutside == false){
	  randomNum = (float)rand()/RAND_MAX;
	  //randomly accept the new position
	  // std::cout << "accept? " << accept  << " random: "<< randomNum << std::endl;	     
	  if(randomNum < accept)
	    {
	      //accept
	      isAccepted = true;
	      //increase density at the new position
	      pDensity.at(newVoxelPos[X_POS],newVoxelPos[Y_POS],newVoxelPos[Z_POS])++;
	      Z_sum[newVoxelPos[Z_POS]]++;
	      //new position becomes old position 
	      (*itr)->SetOldPos(new_pos);
	      if(oldOutside==false){
		//decrease density at old position
		pDensity.at(oldVoxelPos[X_POS],oldVoxelPos[Y_POS],oldVoxelPos[Z_POS])--;
		Z_sum[newVoxelPos[Z_POS]]--;
	      }
	      
	    }
	
	}
      }
  }
  

  //final state:
  std::cout << "making final state" << std::endl;
 for(itr=gammaSet.begin(); itr!=gammaSet.end(); itr++)
      {	
	//get old positions/voxels:
	old_pos[X_POS]=(*itr)->GetOldX();
	old_pos[Y_POS]=(*itr)->GetOldY();
	old_pos[Z_POS]=(*itr)->GetOldZ();

	oldVoxelPos[X_POS] = int((old_pos[X_POS] - xmin)/voxel_dimX);
	oldVoxelPos[Y_POS] = int((old_pos[Y_POS] - ymin)/voxel_dimY);
	oldVoxelPos[Z_POS] = int((old_pos[Z_POS] - zmin)/voxel_dimZ);

	
	if ((oldVoxelPos[X_POS] >= 0)&&(oldVoxelPos[Y_POS] >= 0)&&(oldVoxelPos[Z_POS] >= 0)&&(oldVoxelPos[X_POS] < Xbins) && (oldVoxelPos[Y_POS]<Ybins) && (oldVoxelPos[Z_POS]<Zbins)){
	  
	finalImage.at(oldVoxelPos[X_POS],oldVoxelPos[Y_POS],oldVoxelPos[Z_POS])++;
	}
      }


 std::cout << "finishing final state" << std::endl;


  std::cout << "End of iterations" << std::endl;
  // //write density matrix in file
  
  std::ofstream outFile ("density.img", ios::binary | ios::out);
  int temp;  
   TFile *f1 = new TFile("ntuple.root","recreate");
  //TTree *t1 = new TTree("t1","");

  float zpos,ypos,xpos;
  TNtuple *ntuple = new TNtuple("ntuple","ntuple","xpos:ypos:zpos:temp");

  //  TBranch *b1 = t1->Branch("xitr",&xitr,"xitr/I");
  //TBranch *b2 = t1->Branch("yitr",&yitr,"yitr/I");
  //TBranch *b3 = t1->Branch("zitr",&zitr,"zitr/I");
  //TBranch *b4 = t1->Branch("temp",&temp,"temp/I");

  //test:

  float *voxelImage = new float[Zbins*Ybins*Xbins];

   float sum = 0;
  for (int zitr =0;zitr < Zbins;zitr++){
    sum = 0;
    zpos = (zitr*voxel_dimZ) + zmin;
    for (int yitr =0;yitr < Ybins;yitr++){
      ypos = (yitr*voxel_dimY) + ymin;
      for (int xitr =0;xitr < Xbins;xitr++){
	temp = pDensity.at(xitr,yitr,zitr);
	//temp = finalImage(xitr,yitr,zitr);
	xpos = (xitr*voxel_dimX) + xmin;
	// //voxelImage[zitr*Xbins*Xbins + Ybins*Xbins + xitr] =  temp;
	outFile.write((char*)&temp,sizeof(temp));
	if (temp>0){
    ntuple->Fill(xpos,ypos,zpos,temp);
	  //std::cout <<temp << "\t" <<  xpos << "\t" << ypos << "\t" << zpos << std::endl;
	}
	sum+=temp;
	// if (temp >0){
	//   std::cout <<temp << "\t voxel " <<  xitr << "\t" << yitr << "\t" << zitr << std::endl;
	 
	//   }      
      }
    }
    // std::cout << "zitr = " << zitr << "\t" << sum << "\tzpos\t" << zpos<< std::endl; 
  }

 
 // // t1->Draw("xitr:yitr:zitr:temp");
 ntuple->Write();
 f1->Close();
 outFile.close();
 // /*
 // for (int xitr =0;xitr < Xbins;xitr++)
 //   for (int yitr =0;yitr < Ybins;yitr++)
 //     for (int zitr =0;zitr < Zbins;zitr++)
 //       if(voxelImage[zitr*Xbins*Xbins + Ybins*Xbins + xitr] > 0)
	//  {
	//  }
 // */
	//  //std::cout << "test: " << xitr << "\t" << yitr << "\t" <<zitr <<  "\t"<< voxelImage[zitr*Xbins*Xbins + Ybins*Xbins + xitr] << std::endl;
 // //try record in a different way: 
 // /*
 // FILE * outputFile = fopen("testImage.img","wb");
 // fwrite(voxelImage,sizeof(float),Xbins*Ybins*Zbins,outputFile);
 // fclose(outputFile);
 // */

 //  // // De-Allocate memory to prevent memory leak
 //  // for (int i = 0; i < Zbins; ++i) {
 //  //   for (int j = 0; j < Ybins; ++j)
 //  //     delete [] pDensity[i][j];
 //  //   delete [] pDensity[i];
 //  // }
 //  // delete [] pDensity;

 //  // /* for (int zitr =0;zitr < Zbins;zitr++){
 //  //   std::cout << zitr << "  zsum:  "  << Z_sum[zitr] << std::endl;
 //  // }
 //  // */
  
  for(itr=gammaSet.begin(); itr!=gammaSet.end(); itr++){
    delete *itr;
  }

  gammaSet.clear();
  std::cout << "Total time running: " << float((clock() - t_start))/CLOCKS_PER_SEC << std::endl;
  return 0;
}

