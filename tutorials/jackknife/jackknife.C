#include "TF1.h"
#include "TGraphErrors.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TStyle.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TFile.h"

double catenary( double *x, double *par ){
  double x0 = par[0];
  double y0 = par[1];
  double a  = par[2];
  
  double X = x[0];

  return y0 + (cosh(a*(X-x0))-1.0)/a;
}

void jackknife_fit( const char *inputfilename, const char *outputfilename){
  TFile *fout = new TFile( outputfilename, "RECREATE");

  TF1 *cat_func = new TF1( "cat_func", catenary, -1000.0, 1000.0, 3 );
  double startpar[3] = {-1.83, 1.93, 4.3e-3};

  gStyle->SetOptFit();
  
  cat_func->SetParameters(startpar);

  ifstream infile(inputfilename);

  //read the columns x, y, dy from the file, tab-separated:
  double xtemp, ytemp, eytemp;
  vector<double> x, y, ex, ey;
  
  while( infile >> xtemp >> ytemp >> eytemp ){
    x.push_back(xtemp);
    y.push_back(ytemp);
    ex.push_back(0.0);
    ey.push_back(eytemp);
  }

  TGraphErrors *gdata = new TGraphErrors( x.size(), &(x[0]), &(y[0]), &(ex[0]), &(ey[0]) );
  
  gdata->SetMarkerStyle(20);
  gdata->Draw("AP");

  TFitResultPtr fit_result = gdata->Fit( cat_func, "ES" );
  fit_result->Print("V");

  for(int i=0; i<3; i++){
    startpar[i] = cat_func->GetParameter(i);
  }

  

  TRandom3 num(0);

  TH1D *ha = new TH1D("ha","",250,startpar[2]-5.0*cat_func->GetParError(2), startpar[2]+5.0*cat_func->GetParError(2));
  TH1D *hx0 = new TH1D("hx0","",250,startpar[0]-5.0*cat_func->GetParError(0), startpar[0]+5.0*cat_func->GetParError(0));
  TH1D *hy0 = new TH1D("hy0","",250,startpar[1]-5.0*cat_func->GetParError(1), startpar[1]+5.0*cat_func->GetParError(1));

  TH1D *hchi2 = new TH1D("hchi2","",250,0.0,25.0);
  
  int ntrials = 100000;
  for( int i=0; i<ntrials; i++){
    cat_func->SetParameters(startpar);

    if( (i+1) % 1000 == 0 ) cout << i+1 << endl;
    
    for(int j=0; j<gdata->GetN(); j++){
      gdata->GetPoint(j,xtemp,ytemp);
      eytemp = gdata->GetErrorY(j);
      ytemp = num.Gaus( cat_func->Eval( xtemp ), eytemp );
    
      gdata->SetPoint(j, xtemp, ytemp );
    }

    TFitResultPtr fittemp = gdata->Fit( cat_func, "ESQ" );
    ha->Fill( cat_func->GetParameter(2) );
    hx0->Fill( cat_func->GetParameter(0) );
    hy0->Fill( cat_func->GetParameter(1) );
    hchi2->Fill( fittemp->Chi2() );
    
  }

  gdata->Write("gdata");

  fout->Write();
}
