import sys

from ROOT import *
import ROOT
from array import array
from constants_theory import *
from theory import *

def Crossing(fit, thresh):
    xThresh=fit.GetX(thresh)
    print "Thresh=", thresh, "==== ", (2*xThresh)**2*100

    return (2*xThresh)**2*100


if __name__ == '__main__':
    
    gROOT.ProcessLine(".L ~/CLICdpStyle/rootstyle/CLICdpStyle.C");
    CLICdpStyle()

#    for i in Thresh:
    device_CLIC_100=Device("n", "p", 100e-4, 20, 16, 1160, 10.55, -1.173e+04, 25)
    device_CLIC_100.calculateGraphs()

    device_CLIC_50=Device("n", "p", 50e-4, 20, 16, 1160, 10.55, -1.173e+04, 25)
    device_CLIC_50.calculateGraphs()

    device_TPX_50=Device("n", "p", 50e-4, 20, 16, 1160, 10.55, -1.173e+04, 55)
    device_TPX_50.calculateGraphs()

    gr3_CLIC_100, gr4_CLIC_100=device_CLIC_100.chargeInPixels_3_4(kRed, kBlack, "Theory")
    gr3_CLIC_50, gr4_CLIC_50=device_CLIC_50.chargeInPixels_3_4(kRed, kBlue, "Theory")
    gr3_TPX_50, gr4_TPX_50=device_TPX_50.chargeInPixels_3_4(kRed, kGreen+2, "Theory")


    canv=TCanvas("charge sharing", "charge sharing")
    mg=TMultiGraph("mg", "mg")
    leg=TLegend(0.5, 0.5, 0.8, 0.85)

    mg.Add(gr4_CLIC_50)
    mg.Add(gr4_CLIC_100)
    mg.Add(gr4_TPX_50)
    
    leg.AddEntry(gr4_CLIC_50, "n-in-p, 25 [#mum] pitch, 50 [#mu] thick", "lp")
    leg.AddEntry(gr4_CLIC_100, "n-in-p, 25 [#mum] pitch, 100 [#mu] thick", "lp")
    leg.AddEntry(gr4_TPX_50, "n-in-p, 55 [#mum] pitch, 50 [#mu] thick", "lp")



    mg.Draw("ALP")
    leg.Draw()

    mg.GetXaxis().SetTitle("Hit position X/pitch")
    mg.GetYaxis().SetTitle("Collected charge [e-]")
    mg.GetYaxis().SetTitleOffset(1.26)

    line1=TLine(canv.GetUxmin(),510, 0.55, 510)
    line1.SetLineColor(kGreen+2)
    line1.Draw("same")

    line2=TLine(canv.GetUxmin(), 826, 0.55, 826)
    line2.SetLineColor(kGreen+2)
    line2.SetLineStyle(2)
    line2.Draw("same")
    
    gPad.Update()
    gPad.Modified()

    bla=raw_input()
