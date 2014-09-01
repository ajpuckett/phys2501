#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TROOT.h"
const double g = 9.81;
const double PI = TMath::Pi();

double tau_a(double *x, double *par){
  double d_A = x[0];
  double M_R = par[0];
  double M_A = par[1];
  double M_B = par[2];
  double R = par[3]; 
  double L = par[4];
  double d_B = par[5];
  double M = M_R + M_A + M_B;

  double I_R = M_R/4.0*(pow(L,2) + pow(R,2)/3.0 );

  double a = M_B/M * d_B + (M_B-M_A + M)/M*L/2.0 - M_A/M * d_A;
  double b = -M_B/M * d_B + (M_A-M_B+M)/M * L/2.0 + M_A/M*d_A;
  
  double I_a = I_R + M_A * pow(d_A,2) + M_B*pow(d_B + L,2);
  double I_b = I_R + M_A * pow(d_A+L,2) + M_B*pow(d_B,2);

  double omega_a = sqrt( M * g * a/I_a );
  double tau_a = 2.0*PI/omega_a;

  return tau_a;
}

double tau_b(double *x, double *par){
  double d_A = x[0];
  double M_R = par[0];
  double M_A = par[1];
  double M_B = par[2];
  double R = par[3]; 
  double L = par[4];
  double d_B = par[5];
  double M = M_R + M_A + M_B;

  double I_R = M_R/4.0*(pow(L,2) + pow(R,2)/3.0 );

  double a = M_B/M * d_B + (M_B-M_A + M)/M*L/2.0 - M_A/M * d_A;
  double b = -M_B/M * d_B + (M_A-M_B+M)/M * L/2.0 + M_A/M*d_A;
  
  double I_a = I_R + M_A * pow(d_A,2) + M_B*pow(d_B + L,2);
  double I_b = I_R + M_A * pow(d_A+L,2) + M_B*pow(d_B,2);

  double omega_b = sqrt( M * g * b/I_b );
  double tau_b = 2.0*PI/omega_b;

  return tau_b;
}

void KatersModel(double L=0.987, double R=1.35, double M_A=2.0, double M_B=0.2, double M_R=0.5, double d_B=0.12){
  gStyle->SetOptStat(0);
  gROOT->ProcessLine(".x ~/rootlogon.C");
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetNdivisions(505,"XYZ");
  gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetLabelSize(0.07,"XYZ");
  gStyle->SetTitleSize(0.07, "XYZ");
  TCanvas *c1 = new TCanvas("c1","c1",800,600);

  TF1 *tau_a_func = new TF1("tau_a_func", tau_a, 0.0, (R-L)/2.0, 6);
  TF1 *tau_b_func = new TF1("tau_b_func", tau_b, 0.0, (R-L)/2.0, 6);

  double parameters[6] = {M_R, M_A, M_B, R, L, d_B};
  tau_a_func->SetParameters(parameters);
  tau_b_func->SetParameters(parameters);

  tau_a_func->SetLineStyle(1);
  tau_a_func->SetLineColor(1);
  tau_b_func->SetLineStyle(9);
  tau_b_func->SetLineColor(2);

  TH2D *hframe = new TH2D("hframe","",200,0.0,(R-L)/2.0,200,1.85,2.25);
  hframe->Draw();
  hframe->GetXaxis()->SetTitle("d_{A} (m)");
  hframe->GetYaxis()->SetTitle("#tau (s)" );

  tau_a_func->Draw("same");
  tau_b_func->Draw("same");

  TLegend *leg1 = new TLegend( 0.5, 0.2, 0.95, 0.5);
  leg1->AddEntry( tau_a_func, "#tau_{a}", "l");
  leg1->AddEntry( tau_b_func, "#tau_{b}", "l");
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->Draw();

}
