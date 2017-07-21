Double_t CauchyDens(Double_t *x, Double_t *par)
{
   Double_t pi   = TMath::Pi();
   Double_t mean = par[0];
   Double_t fwhm = par[1];

   Double_t arg = x[0]-mean;
   Double_t top = 0.5*fwhm;
   Double_t bot = pi*(arg*arg+top*top);

   Double_t func = top/bot;
   return func;
}

Double_t CauchyPeak(Double_t *x, Double_t *par)
{
   Double_t height = par[2];
   Double_t func = height*CauchyDens(x,par);
   return func;
}

// Quadratic background function
Double_t background(Double_t *x, Double_t *par) {
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}


// Lorenzian Peak function
Double_t lorentzianPeak(Double_t *x, Double_t *par) {
  return (0.5*par[0]*par[1]/TMath::Pi()) /
    TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2])
   + .25*par[1]*par[1]);
}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
  return background(x,par) + lorentzianPeak(x,&par[3]);
}

void plotImage(){
  Float_t xpos, ypos, zpos, temp;
  TFile *f = new TFile("ntuple.root");
  TNtuple *t1 = (TNtuple*)f->Get("ntuple");
  gStyle->SetOptStat(0);

  
  t1->SetBranchAddress("xpos",&xpos);
  t1->SetBranchAddress("ypos",&ypos);
  t1->SetBranchAddress("zpos",&zpos);
  t1->SetBranchAddress("temp",&temp);
  //t1->Draw("xpos:ypos:zpos:temp","","");


  const Double_t xmin = -200;
  const Double_t xmax = 200;
  const Double_t ymin = -200;
  const Double_t ymax = 200;
  const Int_t nbins = 400;


  TH2D * h2 = new TH2D("","",nbins,xmin,xmax,nbins,ymin,ymax);
   
  Int_t events = 0;
  Int_t sum = 0;
  Int_t nentries = (Int_t)t1->GetEntries();

  for (Int_t i = xmin; i<xmax;i++){
    for (Int_t j = ymin; j<ymax;j++){
      h2->Fill(i,j,1);
    }
  }

  for (Int_t i=0;i<nentries;i++) {
    t1->GetEntry(i);
    sum+=temp;
    
    if (zpos == -10){
      h2->Fill(xpos,ypos,temp);
    }
    if(temp>1)
       events++;
  }
  
  std::cout << "total sum: " << sum << std::endl;
   TCanvas *c1 = new TCanvas("c1","c1",700,700);
  
   h2->Draw("colz");
   //h2->GetZaxis()->SetRangeUser(1,600);
   h2->GetXaxis()->SetTitle("x [mm]");
   h2->GetYaxis()->SetTitle("y [mm]");
   /*
// Creat// Creates a Root function based on function CauchDens above
   TF1 *func = new TF1("cauchy",CauchyDens,-10,10,2);

// Sets initial values and parameter names
   Double_t par[2];
   par[0] = 0.0; par[1] = 10.0;
   func->SetParameters(par);
   func->SetParNames("Mean","FWHM");
// Increase the number of points for calculating the integral of the density
   func->SetNpx(1000);
   Int_t nrPoints = 10000;
// Extract the Cauchy density parameters from the fit
   TF1 *fitfunc = new TF1("fitfunc",CauchyPeak,-100,100,3);
   Double_t fitpar[3];
   fitpar[0] = par[0];
   fitpar[1] = par[1];
   fitpar[2] = Double_t(nrPoints);
   fitfunc->SetParameters(fitpar);
   //hx->Fit(fitfunc);
   */  
   // create a TF1 with the range from 0 to 3 and 6 parameters
   /*TF1 *fitFcn = new TF1("fitFcn",fitFunction,-100,100,6);
   fitFcn->SetNpx(500);
   fitFcn->SetLineWidth(4);
   fitFcn->SetLineColor(kMagenta);

   // first try without starting values for the parameters
   // This defaults to 1 for each param.
   // this results in an ok fit for the polynomial function
   // however the non-linear part (lorenzian) does not
   // respond well.
   fitFcn->SetParameters(1,1,1,1,1,1);
   hx->Fit("fitFcn","0");

   // second try: set start values for some parameters
   fitFcn->SetParameter(4,20); // width
   fitFcn->SetParameter(5,0);   // peak

   TF1* fc = new TF1("fc", "TMath::CauchyDist(x, [0], [1])", -100, 100);
    fc->SetParameters(0, 1);
    //
    hx->Fit("fc");
    //   hx->Fit("fitFcn","V+","ep");
    fc->Draw("same");

   */
   // improve the picture:
   /* TF1 *backFcn = new TF1("backFcn",background,0,3,3);
   backFcn->SetLineColor(kRed);
   TF1 *signalFcn = new TF1("signalFcn",lorentzianPeak,0,3,3);
   signalFcn->SetLineColor(kBlue);
   signalFcn->SetNpx(500);
   Double_t par[6];

   // writes the fit results into the par array
   fitFcn->GetParameters(par);

   backFcn->SetParameters(par);
   backFcn->Draw("same");

   signalFcn->SetParameters(&par[3]);
   signalFcn->Draw("same");
   */
  
}
