void plot_binomial_dist(double p=0.5, int N=100){


  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.045);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadBottomMargin(0.14);

  gStyle->SetTitleFont(62,"XYZ");
  gStyle->SetTextFont(62);
  gStyle->SetLabelFont(62,"XYZ");
  gStyle->SetLegendFont(62);
  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetLabelSize(0.06,"XYZ");

  gStyle->SetOptStat(0);
  //gStyle->SetLineWidth(2);

  gROOT->ForceStyle();

  TH1D *hN5_p = new TH1D("hN5_p","",6,-0.5,5.5);
  TH1D *hN10_p = new TH1D("hN10_p","",11,-0.5,10.5);
  TH1D *hN20_p = new TH1D("hN20_p","",21,-0.5,20.5);
  TH1D *hN50_p = new TH1D("hN50_p","",51,-0.5,50.5);
  TH1D *hN100_p = new TH1D("hN100_p","",101,-0.5,100.5);

  TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
  TCanvas *c2 = new TCanvas("c2","c2",1600,1200);

  for(int i=0; i<=100; i++){
    hN100_p->Fill( i, TMath::Binomial( 100, i ) * pow(p,i)*pow(1.0-p,100-i ) );
    if( i <= 50 ){
      hN50_p->Fill( i, TMath::Binomial( 50, i ) * pow(p,i)*pow(1.0-p,50-i ) );
    }
    if( i <= 20 ){
      hN20_p->Fill( i, TMath::Binomial( 20, i ) * pow(p,i)*pow(1.0-p,20-i) );
    } 
    if( i <= 10 ){
      hN10_p->Fill( i, TMath::Binomial( 10, i ) * pow(p,i)*pow(1.0-p,10-i) );
    }
    if( i <= 5 ){
      hN5_p->Fill( i, TMath::Binomial( 5, i ) * pow(p,i)*pow(1.0-p,5-i) );
    }
  }

  TH1D *hp10_N = new TH1D("hp10_N","",N+1,-0.5,N+0.5);
  TH1D *hp20_N = new TH1D("hp20_N","",N+1,-0.5,N+0.5);
  TH1D *hp50_N = new TH1D("hp50_N","",N+1,-0.5,N+0.5);
  TH1D *hp80_N = new TH1D("hp80_N","",N+1,-0.5,N+0.5);
  TH1D *hp90_N = new TH1D("hp90_N","",N+1,-0.5,N+0.5);

  for(int i=0; i<=N; i++){
    double ptemp = 0.1;
    hp10_N->Fill( i, TMath::Binomial(N,i)*pow(ptemp,i)*pow(1.0-ptemp,N-i) );
    ptemp = 0.2;
    hp20_N->Fill( i, TMath::Binomial(N,i)*pow(ptemp,i)*pow(1.0-ptemp,N-i) );
    ptemp = 0.5;
    hp50_N->Fill( i, TMath::Binomial(N,i)*pow(ptemp,i)*pow(1.0-ptemp,N-i) );
    ptemp = 0.8;
    hp80_N->Fill( i, TMath::Binomial(N,i)*pow(ptemp,i)*pow(1.0-ptemp,N-i) );
    ptemp = 0.9;
    hp90_N->Fill( i, TMath::Binomial(N,i)*pow(ptemp,i)*pow(1.0-ptemp,N-i) );
  }

  c1->cd();

  TH2D *hframe = new TH2D("hframe","",101,-0.5,100.5,100,0.0,1.0);
  hframe->Draw();
  hframe->GetYaxis()->SetTitle("P(k)");
  hframe->GetXaxis()->SetTitle("k");
  
  hp10_N->SetLineWidth(2);
  hp20_N->SetLineWidth(2);
  hp50_N->SetLineWidth(2);
  hp80_N->SetLineWidth(2);
  hp90_N->SetLineWidth(2);

  hN100_p->SetLineWidth(2);
  hN50_p->SetLineWidth(2);
  hN20_p->SetLineWidth(2);
  hN10_p->SetLineWidth(2);
  hN5_p->SetLineWidth(2);

  hN100_p->SetLineStyle(1);
  hN100_p->SetLineColor(1);
  hN100_p->Draw();
  hN100_p->GetYaxis()->SetRangeUser(1.0e-6,1.0);
  hN100_p->GetYaxis()->SetTitle("P(k)");
  hN100_p->GetXaxis()->SetTitle("k");

  hN50_p->SetLineStyle(10);
  hN50_p->SetLineColor(2);
  hN50_p->Draw("same");
  
  hN20_p->SetLineStyle(9);
  hN20_p->SetLineColor(4);
  hN20_p->Draw("same");

  hN10_p->SetLineStyle(8);
  hN10_p->SetLineColor(6);
  hN10_p->Draw("same");
  
  hN5_p->SetLineStyle(6);
  hN5_p->SetLineColor(8);
  hN5_p->Draw("same");

  TLegend *leg = new TLegend(0.5,0.5,0.955, 0.94);

  TString header;

  header.Form("p = %5.3g",p);

  leg->SetHeader( header );
  leg->AddEntry(hN100_p, "N = 100", "l" );
  leg->AddEntry(hN50_p, "N = 50", "l" );
  leg->AddEntry(hN20_p, "N = 20", "l" );
  leg->AddEntry(hN10_p, "N = 10", "l" );
  leg->AddEntry(hN5_p, "N = 5", "l" );
  
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw();

  c2->cd();
  hp10_N->Draw();
  // hp10_N->GetYaxis()->SetRangeUser(1.0e-6,1.0);
  hp10_N->GetYaxis()->SetTitle("P(k)");
  hp10_N->GetXaxis()->SetTitle("k");
  hp10_N->GetYaxis()->SetTitleOffset(1.2);
  hp10_N->SetLineStyle(1);
  hp10_N->SetLineColor(1);

  hp20_N->SetLineStyle(10);
  hp20_N->SetLineColor(2);
  hp20_N->Draw("SAME");
  
  hp50_N->SetLineStyle(9);
  hp50_N->SetLineColor(4);
  hp50_N->Draw("SAME");

  hp80_N->SetLineStyle(8);
  hp80_N->SetLineColor(6);
  hp80_N->Draw("SAME");

  hp90_N->SetLineStyle(6);
  hp90_N->SetLineColor(8);
  hp90_N->Draw("SAME");
  
  TLegend *leg2 = new TLegend(0.4,0.6,0.9,0.9);
  leg2->SetHeader( TString::Format("N = %d",N) );
  leg2->AddEntry(hp10_N,"p = 0.1", "l");
  leg2->AddEntry(hp20_N,"p = 0.2", "l");
  leg2->AddEntry(hp50_N,"p = 0.5", "l");
  leg2->AddEntry(hp80_N,"p = 0.8", "l");
  leg2->AddEntry(hp90_N,"p = 0.9", "l");
  
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->Draw();
  

}
