/*
 * DreamData.cxx
 *
 *  Created on: 29 Aug 2018
 *      Author: bernhardhohlweger
 */

#include "DreamData.h"
#include "TLegend.h"
#include "TStyle.h"
DreamData::DreamData(const char* particlePair)
    : fName(particlePair),
      fCorrelationFunction(nullptr),
      fCorrelationFunctionSimulation(nullptr),
      fSystematics(nullptr),
      fSysError(nullptr),
      fBaseLine(new TF1(Form("%sBaseLine", particlePair), "pol1", 0, 1000)),
      fXMin(0),
      fXMax(0.5),
      fYMin(0),
      fYMax(0.5),
      fXMinLegend(0),
      fXMaxLegend(0.5),
      fYMinLegend(0),
      fYMaxLegend(0.5),
      fUnitConversionData(1),
      fUnitConversionCATS(1),
      fLegendName(),
      fLegendOption(),
      fFemtoModdeled(),
      fFakeGraph() {
  fFillColors = {
    kGray + 1,
    kRed - 10,
    kBlue - 9,
    kGreen - 8,
    kMagenta - 9,
    kOrange - 9,
    kCyan - 3,
    kYellow - 7
  };
  TColor myColor1;
  fColors = {
    kBlack,         //0
    kRed + 1,//1
    kBlue + 2,//2
    kGreen + 3,//3
    kMagenta + 1,//4
    kOrange - 7,//5
    kCyan + 2,//6
    kYellow + 2,//7
    kWhite,//8
    kGreen - 5,//9
    myColor1.GetColor(255,127,0),
    myColor1.GetColor(31,120,180),
    myColor1.GetColor(178,223,138)
  };
  fMarkers = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare, kOpenDiamond,
    kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};
  fBaseLine->SetLineStyle(2);
  fBaseLine->SetLineColor(fFillColors[6]);

}

DreamData::~DreamData() {
  // TODO Auto-generated destructor stub
}

void DreamData::SetSystematics(TF1* parameters, float errorwidth) {
  if (parameters) {
    if (fCorrelationFunction) {
      int nBinsX = fCorrelationFunction->GetNbinsX();
      float minX = fCorrelationFunction->GetXaxis()->GetXmin();
      float maxX = fCorrelationFunction->GetXaxis()->GetXmax();
      fSystematics = new TH1F(Form("Systematics%s", fName),
                              Form("Systematics%s", fName), nBinsX, minX, maxX);
      for (int iBin = 1; iBin < nBinsX; iBin++) {
        const float x = fCorrelationFunction->GetBinCenter(iBin);
        const float y = fCorrelationFunction->GetBinContent(iBin);
//        std::cout << "x = " << x << '\t' << " y = " << y << '\n'
//                  << " rel error "
//                  << parameters->Eval(x / (float) fUnitConversionData) << '\t'
//                  << " Error = "
//                  << y * parameters->Eval(x / (float) fUnitConversionData)
//                  << std::endl;
        fSystematics->SetBinContent(
            iBin, y * parameters->Eval(x / (float) fUnitConversionData));
      }
      fSystematics->SetLineWidth(2.0);
      fSysError = new TGraphErrors();

      for (int i = 0; i < nBinsX; i++) {
        if (fCorrelationFunction->GetBinCenter(i + 1)
            > 0.2 * fUnitConversionData)
          continue;
        fSysError->SetPoint(i, fCorrelationFunction->GetBinCenter(i + 1),
                            fCorrelationFunction->GetBinContent(i + 1));
        fSysError->SetPointError(i, errorwidth,
                                 fSystematics->GetBinContent(i + 1));
      }
      TGraph *grFakeSys = new TGraph();
      SetStyleGraph(grFakeSys, 2, 0);
      grFakeSys->SetFillColor(fFillColors[0]);
      grFakeSys->SetLineColor(fFillColors[0]);
      grFakeSys->SetLineWidth(0);
      fFakeGraph.push_back(grFakeSys);
    } else {
      std::cout << "For " << fName
                << " set the CF before adding the systematics \n";
    }
  } else {
    std::cout << "Paramters input missing for " << fName << std::endl;
  }
  return;
}

void DreamData::FemtoModelFitBands(TGraph *grMedian1, TGraph *grLower,
                                   TGraph *grUpper, int color, int lineStyle,
                                   int lineWidth, int fillStyle,
                                   bool addtoLegend) {
  if (fSystematics) {
    TGraphErrors *grFemtoModel = new TGraphErrors();
    grFemtoModel->SetName(grMedian1->GetName());
    double x, yM1, yLo, yUp;
    int count = 0;
    for (int i = 0; i < grMedian1->GetN(); ++i) {
      grMedian1->GetPoint(i, x, yM1);
      grLower->GetPoint(i, x, yLo);
      grUpper->GetPoint(i, x, yUp);
      std::vector<float> yAll;
      yAll.push_back(yM1);
      yAll.push_back(yLo);
      yAll.push_back(yUp);
      std::sort(yAll.begin(), yAll.end());
      grFemtoModel->SetPoint(count, x / (float) fUnitConversionCATS,
                             (yAll[2] + yAll[0]) / 2.f);
      grFemtoModel->SetPointError(count++, 0,
                                  (yAll[2] + yAll[0]) / 2.f - yAll[0]);
    }
    grFemtoModel->SetLineColor(fColors[color]);
    grFemtoModel->SetFillColor(fColors[color]);
    grFemtoModel->SetLineWidth(lineWidth);
    grFemtoModel->SetLineStyle(lineStyle);
    if (fillStyle > 0)
      grFemtoModel->SetFillStyle(fillStyle);
    fFemtoModdeled.push_back(grFemtoModel);
    if (addtoLegend) {
      TGraph *grFakeModel = new TGraph();
      grFakeModel->SetLineColor(fColors[color]);
      grFakeModel->SetFillColor(fColors[color]);
      grFakeModel->SetLineWidth(lineWidth * 1.8);
      grFakeModel->SetLineStyle(lineStyle);
      if (fillStyle > 0) {
        grFakeModel->SetFillStyle(fillStyle);
      }
      fFakeGraph.push_back(grFakeModel);
    }
  } else {
    std::cout << "Set Systematics first for " << fName << "\n";
  }
  return;
}

void DreamData::SetStyleHisto(TH1 *histo, int marker, int color) {
  histo->GetXaxis()->SetLabelSize(0.045);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.045);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->SetMarkerSize(1.4);
  histo->SetLineWidth(2);
  histo->SetMarkerStyle(fMarkers[marker]);
  histo->SetMarkerColor(fColors[color]);
  histo->SetLineColor(fColors[color]);
}

void DreamData::DrawCorrelationPlot(TCanvas* c) {
  c->cd();
  SetStyleHisto(fCorrelationFunction, 2, 0);
  fCorrelationFunction->GetXaxis()->SetRangeUser(fXMin, fXMax);
  fCorrelationFunction->GetYaxis()->SetRangeUser(fYMin, fYMax);
  fSysError->SetLineColor(kWhite);
  fSysError->Draw("Ap");
  fBaseLine->Draw("same");
  if (fUnitConversionData == 1) {
    fSysError->SetTitle("; #it{k}* (GeV/#it{c}); #it{C}(#it{k}*)");
  } else if (fUnitConversionData) {
    fSysError->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  } else {
    std::cout << "in DreamData::DrawCorrelationPlot we don't know this unit \n";
  }
  fSysError->GetXaxis()->SetRangeUser(fXMin, fXMax);
  fSysError->GetYaxis()->SetRangeUser(fYMin, fYMax);
  TLegend *leg = new TLegend(fXMinLegend, fYMinLegend, fXMaxLegend,
                             fYMaxLegend);
//  TLegend *leg = new TLegend(0.5, 0.55, 0.62, 0.875);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(gStyle->GetTextSize() * 0.75);
  int legendCounter = 1;
//  leg->AddEntry(fCorrelationFunction, fLegendName[0], "pe");
  leg->AddEntry(fFakeGraph[0], fLegendName[0], fLegendOption[0]);
  leg->AddEntry(fBaseLine, "Baseline", "l");
  leg->Draw("same");
  for (auto &it : fFemtoModdeled) {
    std::cout << "leg counter " << legendCounter << std::endl;
    it->Draw("L3 same");
    std::cout << "size " << fFakeGraph.size()<< std::endl;
    if (legendCounter < fFakeGraph.size()) {
      std::cout << " leg counter 2 " << legendCounter << std::endl;
    leg->AddEntry(fFakeGraph[legendCounter], fLegendName[legendCounter],
                  fLegendOption[legendCounter]);
    }
    legendCounter++;
  }
  fSysError->SetFillColorAlpha(kBlack, 0.4);
  fSysError->Draw("2 same");
  fCorrelationFunction->Draw("pe same");
  leg->Draw("same");
}

void DreamData::SetStyleGraph(TGraph *histo, int marker, int color) {
  histo->GetXaxis()->SetLabelSize(0.045);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.045);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->SetMarkerSize(1.4);
  histo->SetLineWidth(2);
  histo->SetMarkerStyle(fMarkers[marker]);
  histo->SetMarkerColor(fColors[color]);
  histo->SetLineColor(fColors[color]);
}
