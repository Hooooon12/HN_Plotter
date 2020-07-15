#include <iostream>

// How to run this code?
// root -b -l -q makePlots_muon_eachYear.C

TString workdir = "/data6/Users/helee/working_HN_Plotter/";
TString SKFlatVersion = "Run2Legacy_v4";
TString skim = "SkimTree_Dilepton";
TString analyzer = "HNtypeI_SR";
TString file_path = "";
vector<TString> year = {"2016", "2017", "2018"};
vector<TString> luminosity = {"35.9", "41.5", "59.7"};
vector<TString> ZGname = {"ZGTo2LG", "ZGToLLG_01J", "ZGToLLG_01J"};
vector<TString> WGname = {"WGToLNuG", "WGToLNuG_01J", "WGToLNuG_01J"};

const int MCNumber = 25;
int maxBinNumber_total = 0, maxBinNumber_temp = 0;
double minRange = 0., maxRange = 0., binContent = 0., binError = 0., binError_Stat = 0., binError_Syst = 0.;
double max_Data = 0., max_Background = 0., max_Hist = 0.;

void FixOverflows(TH1D *hist, int maxBin, int maxBin_total);

void makePlots_muon_eachYear(){

  string histline;
  ifstream in("histList_emu.txt");
  // Line loop
  while(getline(in, histline)){
    std::istringstream is(histline);
    TString IDname, channel, region, variable, txt_channel, PDname, txt_region, output_region, txt_variable;
    int rebin, minBinNumber, maxBinNumber;
    is >> IDname;
    is >> channel;
    is >> region;
    is >> variable;
    is >> rebin;
    is >> minBinNumber;
    is >> maxBinNumber;

    // txt_channel, PDname
    if(channel == "diel"){
      txt_channel = "e^{#pm}e^{#pm}";
      PDname = "DoubleEG";
    }
    if(channel == "dimu"){
      txt_channel = "#mu^{#pm}#mu^{#pm}";
      PDname = "DoubleMuon";
    }
    if(channel == "emu"){
      txt_channel = "e^{#pm}#mu^{#pm}";
      PDname = "MuonEG";
    }

    // txt_region, output_region
    if(region == "Pre") txt_region = "Preselection";
    if(region == "lowCR1") txt_region = "Low Mass CR1";
    if(region == "lowSR1") txt_region = "Low Mass SR1";
    if(region == "lowCR2") txt_region = "Low Mass CR2";
    if(region == "lowSR2") txt_region = "Low Mass SR2";
    if(region == "highCR1") txt_region = "High Mass CR1";
    if(region == "highSR1") txt_region = "High Mass SR1";
    if(region == "highCR2") txt_region = "High Mass CR2";
    if(region == "highSR2") txt_region = "High Mass SR2";
    if(region == "fakeCR1") txt_region = "Non-prompt CR1";
    if(region == "fakeCR2") txt_region = "Non-prompt CR2";
    output_region = region;

    if(variable.Contains("nocut")){ // Preselection + jet requirement for each low/high mass region
       txt_region = "Preselection";
       output_region = "Pre";
    }

    // txt_variable
    if(variable.Contains("Number_Jets")) txt_variable = "Number of AK4 jets";
    if(variable.Contains("Number_FatJets")) txt_variable = "Number of AK8 jets";
    if(variable.Contains("Number_BJets")) txt_variable = "Number of b jets";
    if(variable.Contains("Lep1_Pt")) txt_variable = "p_{T}(l_{1}) (GeV)";
    if(variable.Contains("Lep2_Pt")) txt_variable = "p_{T}(l_{2}) (GeV)";
    if(variable.Contains("MET")) txt_variable = "#slash{E}_{T}^{miss} (GeV)";
    if(variable.Contains("MET2ST")) txt_variable = "(#slash{E}_{T}^{miss})^{2}/S_{T} (GeV)";
    if(variable.Contains("WCand1")) txt_variable = "m(lljj) (GeV)";
    if(variable.Contains("WCand2")) txt_variable = "m(jj) (GeV)";
    if(variable.Contains("llj")) txt_variable = "m(llj) (GeV)";
    if(variable.Contains("l1j")) txt_variable = "m(l_{1}j) (GeV)";
    if(variable.Contains("l2j")) txt_variable = "m(l_{2}j) (GeV)";
    if(variable.Contains("lljj")) txt_variable = "m(lljj) (GeV)";
    if(variable.Contains("l1jj")) txt_variable = "m(l_{1}jj) (GeV)";
    if(variable.Contains("l2jj")) txt_variable = "m(l_{2}jj) (GeV)";
    if(variable.Contains("l1J")) txt_variable = "m(l_{1}J) (GeV)";
    if(variable.Contains("l2J")) txt_variable = "m(l_{2}J) (GeV)";
    if(variable.Contains("FatJet_Pt")) txt_variable = "p_{T}(J) (GeV)";
    if(variable.Contains("FatJet_Mass")) txt_variable = "m(J) (GeV)";

    // Declare variables needed for making plots 
    TFile *f_Data[3], *f_Fake[3], *f_MC[MCNumber][3];
    TH1D *h_Data[3], *h_Fake[3], *h_Temp[3], *h_AllMC[3], *h_MC[MCNumber][3], *h_Error[3], *h_Error_Background1[3], *h_Error_Background2[3], *h_Ratio[3];
    TCanvas *c1;
    TPad *c_up, *c_down;
    THStack *hs;
    TLegend *lg, *lg2;

    // Year loop
    for(int it_y=0; it_y<year.size(); it_y++){
      file_path = SKFlatVersion+"/"+analyzer+"/"+year.at(it_y)+"/";

      // PDname in 2018 : DoubleEG -> EGamma
      if(channel == "diel"){
        if(it_y == 2) PDname = "EGamma";
        else PDname = "DoubleEG";
      }

      //=========================================
      //==== Set input ROOT files
      //=========================================

      // DATA, Fake
      f_Data[it_y]   = new TFile(workdir+file_path+"DATA/"+analyzer+"_"+skim+"_"+PDname+".root");
      f_Fake[it_y]   = new TFile(workdir+file_path+"RunFake__/DATA/"+analyzer+"_"+skim+"_"+PDname+".root");
      // MC : VV, VG
      f_MC[0][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim+"_WZTo3LNu_powheg.root"); 
      f_MC[1][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim+"_ZZTo4L_powheg.root");
      f_MC[2][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim+"_"+ZGname.at(it_y)+".root");
      f_MC[3][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim+"_"+WGname.at(it_y)+".root");
      f_MC[4][it_y]  = new TFile(workdir+file_path+analyzer+"_WWTo2L2Nu_DS.root");
      f_MC[5][it_y]  = new TFile(workdir+file_path+analyzer+"_WpWp_EWK.root");
      f_MC[6][it_y]  = new TFile(workdir+file_path+analyzer+"_WpWp_QCD.root");
      // MC : VVV
      f_MC[7][it_y]  = new TFile(workdir+file_path+analyzer+"_WWW.root");
      f_MC[8][it_y]  = new TFile(workdir+file_path+analyzer+"_WWZ.root");
      f_MC[9][it_y]  = new TFile(workdir+file_path+analyzer+"_WZZ.root");
      f_MC[10][it_y] = new TFile(workdir+file_path+analyzer+"_ZZZ.root");
      // MC : Top
      f_MC[11][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim+"_ttWToLNu.root");
      f_MC[12][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim+"_ttZToLLNuNu.root");
      f_MC[13][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim+"_ttHToNonbb.root");
      f_MC[14][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim+"_TTG.root");
      f_MC[15][it_y] = new TFile(workdir+file_path+analyzer+"_TG.root");
      // MC : Higgs
      f_MC[16][it_y] = new TFile(workdir+file_path+analyzer+"_ggHToZZTo4L.root");
      f_MC[17][it_y] = new TFile(workdir+file_path+analyzer+"_VBF_HToZZTo4L.root");
      f_MC[18][it_y] = new TFile(workdir+file_path+analyzer+"_VHToNonbb.root");
      // MC : ggZZTo4L
      f_MC[19][it_y] = new TFile(workdir+file_path+analyzer+"_ggZZTo2e2mu.root");
      f_MC[20][it_y] = new TFile(workdir+file_path+analyzer+"_ggZZTo2e2tau.root");
      f_MC[21][it_y] = new TFile(workdir+file_path+analyzer+"_ggZZTo2mu2tau.root");
      f_MC[22][it_y] = new TFile(workdir+file_path+analyzer+"_ggZZTo4e.root");
      f_MC[23][it_y] = new TFile(workdir+file_path+analyzer+"_ggZZTo4mu.root");
      f_MC[24][it_y] = new TFile(workdir+file_path+analyzer+"_ggZZTo4tau.root");

      //=========================================
      //==== Get histograms
      //=========================================

      // DATA, Fake
      h_Data[it_y]  = (TH1D*)f_Data[it_y]->Get(channel+"/"+region+"/"+variable+"_"+IDname);
      h_Fake[it_y]  = (TH1D*)f_Fake[it_y]->Get(channel+"/"+region+"/"+variable+"_"+IDname);
      // MC
      for(int it_mc=0; it_mc<MCNumber; it_mc++){
        h_MC[it_mc][it_y] = (TH1D*)f_MC[it_mc][it_y]->Get(channel+"/"+region+"/"+variable+"_"+IDname);
      }

      h_Data[it_y]->SetDirectory(0);
      h_Fake[it_y]->SetDirectory(0);
      for(int it_mc=0; it_mc<MCNumber; it_mc++){
        if(h_MC[it_mc][it_y]) h_MC[it_mc][it_y]->SetDirectory(0);
      }

      f_Data[it_y]->Close();
      f_Fake[it_y]->Close();
      for(int it_mc=0; it_mc<MCNumber; it_mc++){
        f_MC[it_mc][it_y]->Close();
      }

      //=========================================
      //==== Make plots
      //=========================================

      // Make an empty histogram for adding backgrounds
      h_Temp[it_y] = (TH1D*)h_Data[it_y]->Clone();
      maxBinNumber_temp = h_Temp[it_y]->GetNbinsX();
      for(int i=0; i<maxBinNumber_temp+2; i++){
        h_Temp[it_y]->SetBinContent(i, 0.);
        h_Temp[it_y]->SetBinError(i, 0.);
      }

      //==== CANVAS
      c1 = new TCanvas("c1", "", 1000, 1000);
      c1->cd();

      //==== PAD : drawing distribution
      c_up = new TPad("c_up", "", 0, 0.25, 1, 1);
      c_up->SetTopMargin(0.08);
      c_up->SetBottomMargin(0.017);
      c_up->SetLeftMargin(0.14);
      c_up->SetRightMargin(0.04);
      c_up->Draw();
      c_up->cd();

      // Merge all prompt backgrounds (MC)
      h_AllMC[it_y] = (TH1D*)h_Temp[it_y]->Clone();
      for(int it_mc=0; it_mc<MCNumber; it_mc++){
        if(h_MC[it_mc][it_y]) h_AllMC[it_y]->Add(h_MC[it_mc][it_y]);
      }

      // Rebin
      h_Data[it_y]->Rebin(rebin);
      h_Fake[it_y]->Rebin(rebin);
      h_AllMC[it_y]->Rebin(rebin);
      h_Temp[it_y]->Rebin(rebin);

      maxBinNumber_total = h_Data[it_y]->GetNbinsX();  // This is needed for adding overflow bins

      // This is needed for drawing a line in the ratio plot
      minRange = h_Data[it_y]->GetBinLowEdge(minBinNumber);
      maxRange = h_Data[it_y]->GetBinLowEdge(maxBinNumber) + h_Data[it_y]->GetBinWidth(maxBinNumber);

      // Fix overflows
      FixOverflows(h_Data[it_y], maxBinNumber, maxBinNumber_total);
      FixOverflows(h_Fake[it_y], maxBinNumber, maxBinNumber_total);
      FixOverflows(h_AllMC[it_y], maxBinNumber, maxBinNumber_total);

      // Stack & Draw MC
      hs = new THStack("hs", "");
      h_Fake[it_y]->SetLineWidth(0);
      h_Fake[it_y]->SetFillColor(kAzure+8);
      hs->Add(h_Fake[it_y]);
      h_AllMC[it_y]->SetLineWidth(0);
      h_AllMC[it_y]->SetFillColor(kGreen+1);
      hs->Add(h_AllMC[it_y]);
      hs->Draw("hist");
      hs->SetTitle("");
      hs->GetXaxis()->SetLabelSize(0.);
      hs->GetYaxis()->SetLabelSize(0.045);
      hs->GetYaxis()->SetTitle("Events / "+TString::Itoa(rebin, 10)+" GeV");
      hs->GetYaxis()->SetTitleSize(0.075);
      hs->GetYaxis()->SetTitleOffset(0.8);
      hs->GetXaxis()->SetRange(minBinNumber, maxBinNumber);

      // h_Error : A histogram for calcalating the total error of all backgrounds
      h_Error[it_y] = (TH1D*)h_Temp[it_y]->Clone();
      if(h_Fake[it_y]) h_Error[it_y]->Add(h_Fake[it_y]);
      if(h_AllMC[it_y]) h_Error[it_y]->Add(h_AllMC[it_y]); 

      // Add systematic errors
      h_Error_Background1[it_y] = (TH1D*)h_Error[it_y]->Clone();  // Stat. + Syst. // Draw this first in the ratio plot
      h_Error_Background2[it_y] = (TH1D*)h_Error[it_y]->Clone();  // Stat. only
      for(int it_bin = minBinNumber; it_bin < maxBinNumber+1; it_bin++){
        binError_Stat = h_Error_Background1[it_y]->GetBinError(it_bin);
        binError_Syst = h_Fake[it_y]->GetBinContent(it_bin)*0.3;
        binError = sqrt(binError_Stat*binError_Stat + binError_Syst*binError_Syst);
        h_Error[it_y]->SetBinError(it_bin, binError);
      }

      // Draw MC error
      h_Error[it_y]->SetMarkerSize(0);
      h_Error[it_y]->SetLineWidth(0);
      h_Error[it_y]->SetFillStyle(3144);
      h_Error[it_y]->SetFillColor(kBlack);
      h_Error[it_y]->Draw("e2 same");
  
      // Draw Data
      h_Data[it_y]->SetMarkerStyle(20);
      h_Data[it_y]->SetMarkerColor(kBlack);
      h_Data[it_y]->Draw("ep same");

      // Set Min or Max of y axis
      max_Data = h_Data[it_y]->GetBinContent(h_Data[it_y]->GetMaximumBin());
      max_Background  = h_Error[it_y]->GetBinContent(h_Error[it_y]->GetMaximumBin());
      max_Hist = std::max(max_Data, max_Background);
      hs->SetMinimum(0);
      hs->SetMaximum(max_Hist*1.3);

      // Draw the legend
      lg = new TLegend(0.6, 0.45, 0.9, 0.85);
      lg->AddEntry(h_Error[it_y], "Stat. + Syst. Uncertainty", "f");
      lg->AddEntry(h_Data[it_y], "Data", "lep");
      lg->AddEntry(h_AllMC[it_y], "Prompt background", "f");
      lg->AddEntry(h_Fake[it_y], "MisId. Lepton background", "f");
      lg->SetBorderSize(0);
      lg->SetTextSize(0.03);
      lg->SetFillStyle(1001);
      lg->SetShadowColor(0);
      lg->Draw("same");

      // Add text
      TLatex txt;
      txt.SetNDC();
      txt.SetTextSize(0.05);
      txt.SetTextAlign(32);
      txt.SetTextFont(42);
      txt.DrawLatex(.95,.96, luminosity.at(it_y)+" fb^{-1} (13 TeV)");

      TLatex txt2;
      txt2.SetNDC();
      txt2.SetTextSize(0.05);
      txt2.SetTextAlign(12);
      txt2.SetTextFont(62);
      txt2.DrawLatex(.18,.88, txt_channel);

      TLatex txt3;
      txt3.SetNDC();
      txt3.SetTextSize(0.05);
      txt3.SetTextAlign(12);
      txt3.SetTextFont(62);
      txt3.DrawLatex(.18,.82, txt_region);
 
      c1->cd();

      // PAD : drawing ratio
      c_down = new TPad("c_down", "", 0, 0, 1, 0.25);
      c_down->SetTopMargin(0.03);
      c_down->SetBottomMargin(0.35);
      c_down->SetLeftMargin(0.14);
      c_down->SetRightMargin(0.04);
      //c_down->SetGridx();
      //c_down->SetGridy();
      c_down->Draw();
      c_down->cd();

      // Relative error of all backgrounds
      for(int it_bin = minBinNumber; it_bin < maxBinNumber+1; it_bin++){
        binContent = h_Error_Background1[it_y]->GetBinContent(it_bin);
        binError_Stat = h_Error_Background1[it_y]->GetBinError(it_bin);
        binError_Syst = h_Fake[it_y]->GetBinContent(it_bin)*0.3;
        binError  = sqrt(binError_Stat*binError_Stat + binError_Syst*binError_Syst);
        if(binContent != 0.){
          binError = binError/binContent;
          binError_Stat = binError_Stat/binContent;
        }
        else{
          binError = 0.;
          binError_Stat = 0.;
        }
        h_Error_Background1[it_y]->SetBinContent(it_bin, 1.);
        h_Error_Background1[it_y]->SetBinError(it_bin, binError);
        h_Error_Background2[it_y]->SetBinContent(it_bin, 1.);
        h_Error_Background2[it_y]->SetBinError(it_bin, binError_Stat);
      }
      h_Error_Background1[it_y]->SetTitle("");
      h_Error_Background1[it_y]->SetStats(0);
      h_Error_Background1[it_y]->GetXaxis()->SetTitle(txt_variable);
      h_Error_Background1[it_y]->GetYaxis()->SetTitle("#frac{Obs.}{Pred.}");
      h_Error_Background1[it_y]->GetXaxis()->SetRange(minBinNumber, maxBinNumber);
      h_Error_Background1[it_y]->GetYaxis()->SetRangeUser(0., 2.);
      h_Error_Background1[it_y]->GetXaxis()->SetLabelSize(0.12);
      h_Error_Background1[it_y]->GetYaxis()->SetLabelSize(0.08);
      h_Error_Background1[it_y]->GetXaxis()->SetTitleSize(0.16);
      h_Error_Background1[it_y]->GetYaxis()->SetTitleSize(0.14);
      h_Error_Background1[it_y]->GetXaxis()->SetTitleOffset(0.9);
      h_Error_Background1[it_y]->GetYaxis()->SetTitleOffset(0.4);

      h_Error_Background1[it_y]->SetMarkerSize(0);
      h_Error_Background1[it_y]->SetLineWidth(0);
      h_Error_Background1[it_y]->SetFillStyle(1001);
      h_Error_Background1[it_y]->SetFillColor(kGray);
      h_Error_Background1[it_y]->Draw("e2"); 

      h_Error_Background2[it_y]->SetMarkerSize(0);
      h_Error_Background2[it_y]->SetLineWidth(0);
      h_Error_Background2[it_y]->SetFillStyle(1001);
      h_Error_Background2[it_y]->SetFillColor(kCyan);
      h_Error_Background2[it_y]->Draw("e2 same");
 
      // Data/Background ratio
      h_Ratio[it_y] = (TH1D*)h_Data[it_y]->Clone();
      h_Ratio[it_y]->Divide(h_Error[it_y]);
      h_Ratio[it_y]->SetLineColor(1);
      h_Ratio[it_y]->SetMarkerColor(1);
      h_Ratio[it_y]->SetMarkerStyle(20);
      h_Ratio[it_y]->Draw("ep same");

      lg2 = new TLegend(0.75, 0.88, 0.9, 0.95);
      lg2->SetNColumns(2);
      lg2->AddEntry(h_Error_Background2[it_y], "Stat.", "f");
      lg2->AddEntry(h_Error_Background1[it_y], "Stat. + Syst.", "f");
      lg2->SetBorderSize(1);
      lg2->SetTextSize(0.06);
      lg2->SetFillStyle(1001);
      lg2->SetShadowColor(0);
      lg2->Draw("same");

      TLine line(minRange, 1., maxRange, 1.);
      line.SetLineWidth(1);
      line.SetLineColor(2);
      line.Draw();

      c1->cd();

      //=========================================
      //==== Save plots
      //=========================================

      c1->SaveAs("./plots_SR_eachYear/"+IDname+"/"+channel+"/"+output_region+"/"+variable+"_"+year.at(it_y)+".png");

      delete c_up;
      delete c_down;

      c1->Close();

      delete c1;
      delete hs;
      delete lg;
      delete lg2;
      delete f_Data[it_y];
      delete f_Fake[it_y];
      for(int it_mc=0; it_mc<MCNumber; it_mc++){
        delete f_MC[it_mc][it_y];
        delete h_MC[it_mc][it_y];
      }
      delete h_Data[it_y];
      delete h_Fake[it_y];
      delete h_Temp[it_y];
      delete h_AllMC[it_y];
      delete h_Error[it_y];
      delete h_Error_Background1[it_y];
      delete h_Error_Background2[it_y];
      delete h_Ratio[it_y];

    }  // year
  }  // histline  
}


void FixOverflows(TH1D *hist, int maxBin, int maxBin_total){
  double binContent = hist->GetBinContent(maxBin);
  double binError = hist->GetBinError(maxBin)*hist->GetBinError(maxBin);

  for(int i=maxBin+1; i<maxBin_total+2; i++){
    binContent += hist->GetBinContent(i);
    binError += hist->GetBinError(i)*hist->GetBinError(i);
  }

  binError = sqrt(binError);

  hist->SetBinContent(maxBin, binContent);
  hist->SetBinError(maxBin, binError);
}
