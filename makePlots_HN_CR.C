#include <iostream>

// How to run this code?
// root -b -l -q makePlots_SSWW.C

TString workdir = "/data6/Users/jihkim/SKFlatOutput/";
TString SKFlatVersion = "Run2UltraLegacy_v2";
TString skim = "SkimTree_Dilepton";
TString skim2 = "SkimTree_HNMultiLep";
TString analyzer = "Signal";
TString file_path = "";
vector<TString> year = {"2016", "2017", "2018"};
//vector<TString> year = {"2016"};
vector<TString> luminosity = {"36.3", "41.5", "59.8"};
//vector<TString> luminosity = {"36.3"};
//vector<TString> ZGname = {"ZGTo2LG", "ZGToLLG_01J", "ZGToLLG_01J"};
//vector<TString> ZGname = {"ZGToLLG_01J"};
//vector<TString> WGname = {"WGToLNuG", "WGToLNuG_01J", "WGToLNuG_01J"};
//vector<TString> WGname = {"WGToLNuG_01J"};

vector<int> MCNumber = {13, 16, 16};
const int MCNumber_max = 16;
int maxBinNumber_total = 0, maxBinNumber_temp = 0;
double minRange = 0., maxRange = 0., binContent = 0., binError = 0., binError_Stat = 0., binError_Syst = 0.;
double max_Data = 0., max_Background = 0., max_Hist = 0.;

void FixOverflows(TH1D *hist, int maxBin, int maxBin_total);

void makePlots_HN_CR(){

  string histline;
  ifstream in("histList_HN_CR.txt");
  // Line loop
  while(getline(in, histline)){
    std::istringstream is(histline);
    TString this_line = histline;
    if(this_line.Contains("#")||this_line=="") continue;
    TString channel, region, variable, IDname, flag, txt_region, output_region, txt_variable, PDname;
    int rebin, minBinNumber, maxBinNumber;
    is >> channel;
    is >> region;
    is >> variable;
    is >> IDname;
    is >> rebin;
    is >> minBinNumber;
    is >> maxBinNumber;
    is >> flag;

    // txt_region, output_region
    if(channel.Contains("dimu")){
      PDname = "DoubleMuon";
    }
    else if(channel.Contains("diel")){
      PDname = "DoubleEG";
    }
    txt_region = region;
    output_region = region;

    // txt_variable
    if(variable.Contains("Lep1_Pt")) txt_variable = "p_{T}(l_{1}) (GeV)";
    if(variable.Contains("Lep2_Pt")) txt_variable = "p_{T}(l_{2}) (GeV)";
    if(variable.Contains("Lep3_Pt")) txt_variable = "p_{T}(l_{3}) (GeV)";
    if(variable.Contains("Lep1_Eta")) txt_variable = "#eta(l_{1})";
    if(variable.Contains("Lep2_Eta")) txt_variable = "#eta(l_{2})";
    if(variable.Contains("Mt")) txt_variable = "M_{T}(l,#slash{E}_{T}^{miss}) (GeV)";
    if(variable.Contains("MET")) txt_variable = "#slash{E}_{T}^{miss} (GeV)";
    if(variable.Contains("MET2ST")) txt_variable = "(#slash{E}_{T}^{miss})^{2}/S_{T} (GeV)";
    if(variable.Contains("ZCand")) txt_variable = "m(ll) (GeV)";
    if(variable.Contains("TriLep")) txt_variable = "m(lll) (GeV)";
    if(variable.Contains("l1jj_Mass")) txt_variable = "m(l_{1}jj) (GeV)";
    if(variable.Contains("l2jj_Mass")) txt_variable = "m(l_{2}jj) (GeV)";
    if(variable.Contains("l1J_Mass")) txt_variable = "m(l_{1}J) (GeV)";
    if(variable.Contains("l2J_Mass")) txt_variable = "m(l_{2}J) (GeV)";
    if(variable.Contains("Fatjet_Mass")) txt_variable = "m(J) (GeV)";

    // Declare variables needed for making plots 
    TFile *f_Data[3], *f_Fake[3], *f_MC[MCNumber_max][3];
    TH1D *h_Data[3], *h_Fake[3], *h_Temp[3], *h_Bundle[3][3], *h_MC[MCNumber_max][3], *h_Error[3], *h_Error_Background1[3], *h_Error_Background2[3], *h_Ratio[3];
    TCanvas *c1;
    TPad *c_up, *c_down;
    THStack *hs;
    TLegend *lg, *lg2;

    // Year loop
    for(int it_y=0; it_y<year.size(); it_y++){
      file_path = SKFlatVersion+"/"+analyzer+"/"+year.at(it_y)+"/";

      // PDname in 2018 : DoubleEG -> EGamma
      if(channel.Contains("diel")){
        if(year[it_y]=="2018") PDname = "EGamma";
        else PDname = "DoubleEG";
      }

      //=========================================
      //==== Set input ROOT files
      //=========================================

      // DATA, Fake
      f_Data[it_y]   = new TFile(workdir+file_path+"DATA/"+analyzer+"_"+skim+"_"+PDname+".root");
      if(flag=="FR_ex") f_Fake[it_y] = new TFile(workdir+file_path+"RunFake__FR_ex__/DATA/"+analyzer+"_"+skim+"_"+PDname+".root"); //JH
      else f_Fake[it_y] = new TFile(workdir+file_path+"RunFake__/DATA/"+analyzer+"_"+skim+"_"+PDname+".root"); //JH
      //MC : WW
      if(year[it_y]=="2016"){
        f_MC[0][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WpWpJJ_EWK.root"); //JH
        f_MC[1][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WpWpJJ_QCD.root");
      }
      else{
        f_MC[0][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WpWpJJ_EWK.root"); //JH
        f_MC[1][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WpWpJJ_QCD.root");
        f_MC[2][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WWTo2L2Nu_DS.root");
      }
      // MC : top
      //f_MC[2][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim2+"_ttHToNonbb.root");
      if(year[it_y]=="2016"){
        f_MC[2][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim2+"_ttWToLNu.root");
        f_MC[3][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim2+"_ttZToLLNuNu.root");
        f_MC[4][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim2+"_SingleTop_tW_top_NoFullyHad.root");
        f_MC[5][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim2+"_SingleTop_tW_antitop_NoFullyHad.root");
      }
      else{
        f_MC[3][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim2+"_ttWToLNu.root");
        f_MC[4][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim2+"_ttZToLLNuNu.root");
        f_MC[5][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim2+"_SingleTop_tW_top_NoFullyHad.root");
        f_MC[6][it_y]  = new TFile(workdir+file_path+analyzer+"_"+skim2+"_SingleTop_tW_antitop_NoFullyHad.root");
      }
      // MC : VVV
      if(year[it_y]=="2016"){
        f_MC[6][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WWW.root");
        f_MC[7][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WWZ.root");
        f_MC[8][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WZZ.root");
        f_MC[9][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_ZZZ.root");
      }
      else{
        f_MC[7][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WWW.root");
        f_MC[8][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WWZ.root");
        f_MC[9][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WZZ.root");
        f_MC[10][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_ZZZ.root");
      }
      // MC : X+G
      if(year[it_y]=="2016"){
        f_MC[10][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_TTG.root");
      }
      else{
        f_MC[11][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_TTG.root");
        f_MC[12][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_ZGTo2LG_01J.root");
        f_MC[13][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WGToLNuG.root");
      }
      //f_MC[13][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_TG.root"); //JH : contains only negative events
      // MC : VV
      if(year[it_y]=="2016"){
        //f_MC[11][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WZTo3LNu_mll0p1_powheg.root");
        f_MC[11][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WZTo3LNu_mllmin4p0_powheg.root");
        f_MC[12][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_ZZTo4L_m_1toInf_powheg.root");
      }
      else{
        //f_MC[11][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WZTo3LNu_mll0p1_powheg.root");
        f_MC[14][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_WZTo3LNu_mllmin4p0_powheg.root");
        f_MC[15][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_ZZTo4L_powheg.root");
        //f_MC[12][it_y] = new TFile(workdir+file_path+analyzer+"_"+skim2+"_ZZTo4L_m_1toInf_powheg.root");
      }

      //=========================================
      //==== Get histograms
      //=========================================

      // DATA, Fake
      h_Data[it_y]  = (TH1D*)f_Data[it_y]->Get(channel+"/"+region+"/"+variable+"_"+IDname);
      h_Fake[it_y]  = (TH1D*)f_Fake[it_y]->Get(channel+"/"+region+"/"+variable+"_"+IDname);
      // MC
      for(int it_mc=0; it_mc<MCNumber[it_y]; it_mc++){
        h_MC[it_mc][it_y] = (TH1D*)f_MC[it_mc][it_y]->Get(channel+"/"+region+"/"+variable+"_"+IDname);
      }

      h_Data[it_y]->SetDirectory(0);
      h_Fake[it_y]->SetDirectory(0);
      for(int it_mc=0; it_mc<MCNumber[it_y]; it_mc++){
        if(h_MC[it_mc][it_y]) h_MC[it_mc][it_y]->SetDirectory(0);
      }

      f_Data[it_y]->Close();
      f_Fake[it_y]->Close();
      for(int it_mc=0; it_mc<MCNumber[it_y]; it_mc++){
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
      //c_up->SetLogy(); //JH
      c_up->Draw();
      c_up->cd();

      // Merge backgrounds
      h_Bundle[0][it_y] = (TH1D*)h_Temp[it_y]->Clone(); //JH : WW
      h_Bundle[1][it_y] = (TH1D*)h_Temp[it_y]->Clone(); //JH : top
      h_Bundle[2][it_y] = (TH1D*)h_Temp[it_y]->Clone(); //JH : VVV
      h_Bundle[3][it_y] = (TH1D*)h_Temp[it_y]->Clone(); //JH : X+G
      h_Bundle[4][it_y] = (TH1D*)h_Temp[it_y]->Clone(); //JH : VV
      if(year[it_y]=="2016"){
        for(int it_mc=0; it_mc<=1; it_mc++){
          if(h_MC[it_mc][it_y]) h_Bundle[0][it_y]->Add(h_MC[it_mc][it_y]);
        }
        for(int it_mc=2; it_mc<=5; it_mc++){
          if(h_MC[it_mc][it_y]) h_Bundle[1][it_y]->Add(h_MC[it_mc][it_y]);
        }
        for(int it_mc=6; it_mc<=9; it_mc++){
          if(h_MC[it_mc][it_y]) h_Bundle[2][it_y]->Add(h_MC[it_mc][it_y]);
        }
        for(int it_mc=10; it_mc<=10; it_mc++){
          if(h_MC[it_mc][it_y]) h_Bundle[3][it_y]->Add(h_MC[it_mc][it_y]);
        }
        for(int it_mc=11; it_mc<MCNumber[it_y]; it_mc++){
          if(h_MC[it_mc][it_y]) h_Bundle[4][it_y]->Add(h_MC[it_mc][it_y]);
        }
      }
      else{
        for(int it_mc=0; it_mc<=2; it_mc++){
          if(h_MC[it_mc][it_y]) h_Bundle[0][it_y]->Add(h_MC[it_mc][it_y]);
        }
        for(int it_mc=3; it_mc<=6; it_mc++){
          if(h_MC[it_mc][it_y]) h_Bundle[1][it_y]->Add(h_MC[it_mc][it_y]);
        }
        for(int it_mc=7; it_mc<=10; it_mc++){
          if(h_MC[it_mc][it_y]) h_Bundle[2][it_y]->Add(h_MC[it_mc][it_y]);
        }
        for(int it_mc=11; it_mc<=13; it_mc++){
          if(h_MC[it_mc][it_y]) h_Bundle[3][it_y]->Add(h_MC[it_mc][it_y]);
        }
        for(int it_mc=14; it_mc<MCNumber[it_y]; it_mc++){
          if(h_MC[it_mc][it_y]) h_Bundle[4][it_y]->Add(h_MC[it_mc][it_y]);
        }
      }

      // Rebin
      h_Data[it_y]->Rebin(rebin);
      h_Fake[it_y]->Rebin(rebin);
      //h_MC[0][it_y]->Rebin(rebin);
      //h_MC[1][it_y]->Rebin(rebin);
      //h_MC[2][it_y]->Rebin(rebin);
      h_Bundle[0][it_y]->Rebin(rebin);
      h_Bundle[1][it_y]->Rebin(rebin);
      h_Bundle[2][it_y]->Rebin(rebin);
      h_Bundle[3][it_y]->Rebin(rebin);
      h_Bundle[4][it_y]->Rebin(rebin);
      h_Temp[it_y]->Rebin(rebin);

      maxBinNumber_total = h_Data[it_y]->GetNbinsX();  // This is needed for adding overflow bins

      // This is needed for drawing a line in the ratio plot
      minRange = h_Data[it_y]->GetBinLowEdge(minBinNumber);
      maxRange = h_Data[it_y]->GetBinLowEdge(maxBinNumber) + h_Data[it_y]->GetBinWidth(maxBinNumber);

      // Fix overflows
      FixOverflows(h_Data[it_y], maxBinNumber, maxBinNumber_total);
      FixOverflows(h_Fake[it_y], maxBinNumber, maxBinNumber_total);
      //FixOverflows(h_MC[0][it_y], maxBinNumber, maxBinNumber_total);
      //FixOverflows(h_MC[1][it_y], maxBinNumber, maxBinNumber_total);
      //FixOverflows(h_MC[2][it_y], maxBinNumber, maxBinNumber_total); //JH
      FixOverflows(h_Bundle[0][it_y], maxBinNumber, maxBinNumber_total);
      FixOverflows(h_Bundle[1][it_y], maxBinNumber, maxBinNumber_total);
      FixOverflows(h_Bundle[2][it_y], maxBinNumber, maxBinNumber_total); //JH
      FixOverflows(h_Bundle[3][it_y], maxBinNumber, maxBinNumber_total); //JH
      FixOverflows(h_Bundle[4][it_y], maxBinNumber, maxBinNumber_total); //JH

      // Stack & Draw MC
      hs = new THStack("hs", "");
      //h_MC[0][it_y]->SetLineWidth(0);
      //h_MC[0][it_y]->SetFillColor(kGreen);
      //hs->Add(h_MC[0][it_y]);
      //h_MC[1][it_y]->SetLineWidth(0);
      //h_MC[1][it_y]->SetFillColor(kRed-1);
      //hs->Add(h_MC[1][it_y]);
      //h_MC[2][it_y]->SetLineWidth(0);
      //h_MC[2][it_y]->SetFillColor(kMagenta);
      //hs->Add(h_MC[2][it_y]);
      h_Fake[it_y]->SetLineWidth(0);
      h_Fake[it_y]->SetFillColor(kAzure+8);
      hs->Add(h_Fake[it_y]);
      h_Bundle[4][it_y]->SetLineWidth(0);
      h_Bundle[4][it_y]->SetFillColor(kSpring+10);
      hs->Add(h_Bundle[4][it_y]);
      h_Bundle[3][it_y]->SetLineWidth(0);
      h_Bundle[3][it_y]->SetFillColor(kGreen+3);
      hs->Add(h_Bundle[3][it_y]);
      h_Bundle[2][it_y]->SetLineWidth(0);
      h_Bundle[2][it_y]->SetFillColor(kOrange);
      hs->Add(h_Bundle[2][it_y]);
      h_Bundle[1][it_y]->SetLineWidth(0);
      h_Bundle[1][it_y]->SetFillColor(kRed);
      hs->Add(h_Bundle[1][it_y]);
      h_Bundle[0][it_y]->SetLineWidth(0);
      h_Bundle[0][it_y]->SetFillColor(kGreen+1);
      hs->Add(h_Bundle[0][it_y]);
      hs->Draw("hist");
      hs->SetTitle("");
      hs->GetXaxis()->SetLabelSize(0.);
      hs->GetYaxis()->SetLabelSize(0.045);
      if(variable.Contains("Eta")) hs->GetYaxis()->SetTitle("Events");
      else hs->GetYaxis()->SetTitle("Events / "+TString::Itoa(rebin, 10)+" GeV");
      hs->GetYaxis()->SetTitleSize(0.075);
      hs->GetYaxis()->SetTitleOffset(0.8);
      hs->GetXaxis()->SetRange(minBinNumber, maxBinNumber);

      // h_Error : A histogram for calcalating the total error of all backgrounds
      h_Error[it_y] = (TH1D*)h_Temp[it_y]->Clone();
      if(h_Fake[it_y]) h_Error[it_y]->Add(h_Fake[it_y]);
      //if(h_MC[0][it_y]) h_Error[it_y]->Add(h_MC[0][it_y]);
      //if(h_MC[1][it_y]) h_Error[it_y]->Add(h_MC[1][it_y]);
      //if(h_MC[2][it_y]) h_Error[it_y]->Add(h_MC[2][it_y]);
      if(h_Bundle[0][it_y]) h_Error[it_y]->Add(h_Bundle[0][it_y]);
      if(h_Bundle[1][it_y]) h_Error[it_y]->Add(h_Bundle[1][it_y]);
      if(h_Bundle[2][it_y]) h_Error[it_y]->Add(h_Bundle[2][it_y]);
      if(h_Bundle[3][it_y]) h_Error[it_y]->Add(h_Bundle[3][it_y]);
      if(h_Bundle[4][it_y]) h_Error[it_y]->Add(h_Bundle[4][it_y]);

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
      max_Background = h_Error[it_y]->GetBinContent(h_Error[it_y]->GetMaximumBin());
      max_Hist = std::max(max_Data, max_Background);
      hs->SetMinimum(0);
      //hs->SetMaximum(max_Hist*10); //JH : use this when using SetLogy
      hs->SetMaximum(max_Hist+30); //JH

      // Draw the legend
      if(variable.Contains("Eta")) lg = new TLegend(0.65, 0.65, 0.9, 0.9);
      else lg = new TLegend(0.6, 0.45, 0.9, 0.85);
      lg->AddEntry(h_Error[it_y], "Stat. + Syst. Uncertainty", "f");
      lg->AddEntry(h_Data[it_y], "Data", "lep");
      lg->AddEntry(h_Bundle[0][it_y], "WW", "f");
      lg->AddEntry(h_Bundle[1][it_y], "top", "f");
      lg->AddEntry(h_Bundle[2][it_y], "VVV", "f");
      lg->AddEntry(h_Bundle[3][it_y], "ttG", "f");
      lg->AddEntry(h_Bundle[4][it_y], "VV", "f");
      lg->AddEntry(h_Fake[it_y], "MisId. Lepton background", "f");
      //lg->AddEntry(h_MC[2][it_y], "Z#gamma", "f");
      //lg->AddEntry(h_MC[1][it_y], "ZZ", "f");
      //lg->AddEntry(h_MC[0][it_y], "WZ", "f");
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

      /*TLatex txt2;
      txt2.SetNDC();
      txt2.SetTextSize(0.05);
      txt2.SetTextAlign(12);
      txt2.SetTextFont(62);
      txt2.DrawLatex(.18,.88, txt_channel);*/

      TLatex txt3;
      txt3.SetNDC();
      txt3.SetTextSize(0.05);
      txt3.SetTextAlign(12);
      txt3.SetTextFont(62);
      txt3.DrawLatex(.18,.88, txt_region);
 
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
      //h_Error_Background1[it_y]->GetYaxis()->SetRangeUser(0.7, 1.3); //JH
      h_Error_Background1[it_y]->GetYaxis()->SetRangeUser(0, 2);
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
      gSystem->Exec("mkdir -p plots_HN_CR/"+year.at(it_y)+"/plots_HN_"+region);
      if(flag=="FR_ex") c1->SaveAs("./plots_HN_CR/"+year.at(it_y)+"/plots_HN_"+region+"/"+variable+"_"+IDname+"_"+year.at(it_y)+"_FRex.png");
      else c1->SaveAs("./plots_HN_CR/"+year.at(it_y)+"/plots_HN_"+region+"/"+variable+"_"+IDname+"_"+year.at(it_y)+".png");

      delete c_up;
      delete c_down;

      c1->Close();

      delete c1;
      delete hs;
      delete lg;
      delete lg2;
      delete f_Data[it_y];
      delete f_Fake[it_y];
      for(int it_mc=0; it_mc<MCNumber[it_y]; it_mc++){
        delete f_MC[it_mc][it_y];
        delete h_MC[it_mc][it_y];
      }
      delete h_Data[it_y];
      delete h_Fake[it_y];
      delete h_Temp[it_y];
      delete h_Bundle[0][it_y];
      delete h_Bundle[1][it_y];
      delete h_Bundle[2][it_y]; //JH
      delete h_Bundle[3][it_y]; //JH
      delete h_Bundle[4][it_y]; //JH
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
