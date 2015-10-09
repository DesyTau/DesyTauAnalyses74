#include "HttStylesNew.cc"
#include "HtoH.h"

void PlotDimuons(TString histName = "massSelH",
		 TString xtitle = "dimuon mass [GeV]",
		 TString ytitle = "Events",
		 float xLower = 41,
		 float xUpper = 200,
		 bool logY = true) {

  SetStyle();

  TFile * file = new TFile("SingleMu_Run2015B.root");
  
  TString samples[7] = {"ZZ_PY8_50ns",
			"WZ_PY8_50ns",
			"WW_PY8_50ns",
			"WJetsToLNu_aMCatNLO_50ns",
			"TTJets_aMCatNLO_50ns",
			"DYJetsToLL_M-10_50_aMCatNLO_50ns",
			"DYJetsToLL_M-50_aMCatNLO_50ns"
  };

  float xsec[7] = {10.3,22.8,63.2,61526,832,18610,6025};

  float lumi = 40;

  TH1F * histDataOld = (TH1F*)file->Get(histName);

  int nBins = histDataOld->GetNbinsX();
  float xMin = histDataOld->GetBinLowEdge(1);
  float xMax = histDataOld->GetBinLowEdge(nBins+1);

  std::cout << std::endl;
  std::cout << "Histogram " << histName << " : " << "nbins = " << nBins
	    << " , min = " << xMin
	    << " , max = " << xMax << std::endl;
  std::cout << std::endl;

  float bins[100];
  int nBinsNew = nBins;
  std::cout << "New number of bins : ";
  std::cin >> nBinsNew;

  if (nBins % nBinsNew >0) { 
    std::cout << "new number of bins = " << nBinsNew 
	      << "  not multiple of " << nBins << std::endl;
    return;
  }
  float binWidth = (xMax-xMin)/float(nBinsNew);
  for (int iB=0; iB<=nBinsNew; ++iB)
    bins[iB] = xMin + float(iB)*binWidth;

  TH1F * histData = TH1toTH1(histDataOld,nBinsNew,bins,true,"_Data_new");

  TH1F * ewkHist = new TH1F("ewkHist","",nBinsNew,bins);
  TH1F * ttHist  = new TH1F("ttHist","",nBinsNew,bins);
  TH1F * qcdHist = new TH1F("qcdHist","",nBinsNew,bins);
  TH1F * zHist = new TH1F("zHist","",nBinsNew,bins);

  int nSamples = 7;

  //  return;

  for (int iS=0; iS<nSamples; ++iS) {
    //    std::cout << "Sample = " << iS << std::endl;
    TFile * fileMC = new TFile(samples[iS]+".root");
    TH1F * histOld = (TH1F*)fileMC->Get(histName);
    TH1F * hist = TH1toTH1(histOld,nBinsNew,bins,true,"_new_"+samples[iS]);
    TH1F * eventCount = (TH1F*)fileMC->Get("histWeightsH");
    float nGen = eventCount->GetSumOfWeights();
    float norm = xsec[iS]*lumi/nGen;
    TH1F * tempHist = ewkHist;
    if (iS==3)
      tempHist = qcdHist;
    if (iS==4)
      tempHist = ttHist;
    if (iS==5 || iS==6)
      tempHist = zHist;
    tempHist->Add(tempHist,hist,1.,norm);

  }

  //  float dataEvents = 0;
  //  float ttEvents = 0;

  qcdHist->Add(qcdHist,ewkHist);
  ttHist->Add(ttHist,qcdHist);
  zHist->Add(zHist,ttHist);

  TH1F * bkgdErr = (TH1F*)qcdHist->Clone("bkgdErr");
  bkgdErr->SetFillStyle(3013);
  bkgdErr->SetFillColor(1);
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);  
  
  for (int iB=1; iB<=nBinsNew; ++iB) {
    ewkHist->SetBinError(iB,0);
    ttHist->SetBinError(iB,0);
    qcdHist->SetBinError(iB,0);
    zHist->SetBinError(iB,0);
  }

  InitData(histData);
  InitHist(qcdHist,"","",kMagenta,1001);
  InitHist(ttHist,"","",kCyan,1001);
  InitHist(ewkHist,"","",kBlue-4,1001);
  InitHist(zHist,"","",kYellow,1001);
  histData->GetXaxis()->SetTitle(xtitle);
  histData->GetYaxis()->SetTitle(ytitle);
  histData->GetYaxis()->SetTitleOffset(1.3);
  histData->GetYaxis()->SetTitleSize(0.06);
  histData->GetXaxis()->SetRangeUser(xLower,xUpper);
  float yUpper = histData->GetMaximum();
  if (logY)
    histData->GetYaxis()->SetRangeUser(0.5,2*yUpper);
  else
    histData->GetYaxis()->SetRangeUser(0,1.2*yUpper);

  histData->SetMarkerSize(1.5);
  histData->GetXaxis()->SetLabelSize(0);
  histData->GetYaxis()->SetLabelSize(0.06);

  //  nData = histData->GetSum();
  //  float nMC   = ttHist->GetSum();
  //  float eData = TMath::Sqrt(nData);


  TCanvas * canv1 = MakeCanvas("canv1", "", 700, 800);
  TPad *upper = new TPad("upper", "pad",0,0.31,1,1);
  upper->Draw();
  upper->cd();
  upper->SetFillColor(0);
  upper->SetBorderMode(0);
  upper->SetBorderSize(10);
  upper->SetTickx(1);
  upper->SetTicky(1);
  upper->SetLeftMargin(0.17);
  upper->SetRightMargin(0.05);
  upper->SetBottomMargin(0.02);
  upper->SetFrameFillStyle(0);
  upper->SetFrameLineStyle(0);
  upper->SetFrameLineWidth(2);
  upper->SetFrameBorderMode(0);
  upper->SetFrameBorderSize(10);
  upper->SetFrameFillStyle(0);
  upper->SetFrameLineStyle(0);
  upper->SetFrameLineWidth(2);
  upper->SetFrameBorderMode(0);
  upper->SetFrameBorderSize(10);

  histData->Draw("e1");
  zHist->Draw("sameh");
  ttHist->Draw("sameh");
  qcdHist->Draw("sameh");
  ewkHist->Draw("sameh");
  histData->Draw("e1same");

  float chi2 = 0;
  for (int iB=1; iB<=nBinsNew; ++iB) {
    float xData = histData->GetBinContent(iB);
    float xMC = zHist->GetBinContent(iB);
    if (xMC>1e-1) {
      float diff2 = (xData-xMC)*(xData-xMC);
      chi2 += diff2/xMC;
    }
  }
  std::cout << std::endl;
  std::cout << "Chi2 = " << chi2 << std::endl;
  std::cout << std::endl;

  TLegend * leg = new TLegend(0.65,0.6,0.9,0.88);
  SetLegendStyle(leg);
  leg->SetTextSize(0.05);
  leg->AddEntry(histData,"Data","lp");
  leg->AddEntry(ewkHist,"dibosons","f");
  leg->AddEntry(qcdHist,"W+Jets","f");
  leg->AddEntry(ttHist,"t#bar{t}","f");
  leg->AddEntry(zHist,"Z#rightarrow#mu#mu","f");
  leg->Draw();

  TLatex * cms = new TLatex(0.25,0.94,"CMS Preliminary   L = 41 pb^{-1} at #sqrt{s} = 13 TeV");

  cms->SetNDC();
  cms->SetTextSize(0.05);
  cms->Draw();

  if (logY) upper->SetLogy(true);
    
  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv1->cd();

  TH1F * ratioH = (TH1F*)histData->Clone("ratioH");
  ratioH->SetMarkerColor(1);
  ratioH->SetMarkerStyle(20);
  ratioH->SetMarkerSize(1.5);
  ratioH->SetLineColor(1);
  ratioH->GetYaxis()->SetRangeUser(0,2.2);
  ratioH->GetYaxis()->SetNdivisions(505);
  ratioH->GetXaxis()->SetLabelFont(42);
  ratioH->GetXaxis()->SetLabelOffset(0.04);
  ratioH->GetXaxis()->SetLabelSize(0.14);
  ratioH->GetXaxis()->SetTitleSize(0.13);
  ratioH->GetXaxis()->SetTitleOffset(1.2);
  ratioH->GetYaxis()->SetTitle("obs/exp");
  ratioH->GetYaxis()->SetLabelFont(42);
  ratioH->GetYaxis()->SetLabelOffset(0.015);
  ratioH->GetYaxis()->SetLabelSize(0.13);
  ratioH->GetYaxis()->SetTitleSize(0.14);
  ratioH->GetYaxis()->SetTitleOffset(0.5);
  ratioH->GetXaxis()->SetTickLength(0.07);
  ratioH->GetYaxis()->SetTickLength(0.04);
  ratioH->GetYaxis()->SetLabelOffset(0.01);

  for (int iB=1; iB<=nBinsNew; ++iB) {
    float x1 = histData->GetBinContent(iB);
    float x2 = zHist->GetBinContent(iB);
    if (x1>0&&x2>0) {
      float e1 = histData->GetBinError(iB);
      float ratio = x1/x2;
      float eratio = e1/x2;
      ratioH->SetBinContent(iB,ratio);
      ratioH->SetBinError(iB,eratio);
    }
    else {
      ratioH->SetBinContent(iB,1000);
    }
  }


  // ------------>Primitives in pad: lower
  lower = new TPad("lower", "pad",0,0,1,0.30);
  lower->Draw();
  lower->cd();
  lower->SetFillColor(0);
  lower->SetBorderMode(0);
  lower->SetBorderSize(10);
  lower->SetGridy();
  lower->SetTickx(1);
  lower->SetTicky(1);
  lower->SetLeftMargin(0.17);
  lower->SetRightMargin(0.05);
  lower->SetTopMargin(0.026);
  lower->SetBottomMargin(0.35);
  lower->SetFrameFillStyle(0);
  lower->SetFrameLineStyle(0);
  lower->SetFrameLineWidth(2);
  lower->SetFrameBorderMode(0);
  lower->SetFrameBorderSize(10);
  lower->SetFrameFillStyle(0);
  lower->SetFrameLineStyle(0);
  lower->SetFrameLineWidth(2);
  lower->SetFrameBorderMode(0);
  lower->SetFrameBorderSize(10);

  ratioH->Draw("e1");

  lower->Modified();
  lower->RedrawAxis();
  canv1->cd();
  canv1->Modified();
  canv1->cd();
  canv1->SetSelected(canv1);

  canv1->Print(histName+".png");
  canv1->Print(histName+".pdf","Portrait pdf");


}
