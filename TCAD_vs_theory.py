import sys

from ROOT import *
import ROOT
from array import array
from constants_theory import *
from theory import *

def returnTCAD(filename, color):
    f = open(filename)
    lines = [line.strip() for line in open(filename)]
    f.close()

    x=[]
    y=[]

    for i in range(1, len(lines)-1):
        vals=lines[i].split(",")

        x.append((float(vals[0])))
        y.append((float(vals[1])))
        
    graph=TGraph(len(x), array('f', x), array('f', y)) #[V/cm]
    graph.SetMarkerColor(color)
    graph.SetMarkerStyle(7)
    graph.SetLineColor(color)

    return graph

def TCAD_chargeInPixel3_4(colorPixel3, colorPixel4, bias):
    # TCAD
#    TCADpositions=[0, 0.25, 0.3, 0.35, 0.375, 0.4, 0.45, 0.475, 0.5]
    TCADpositions=[0, 0.25, 0.3, 0.35, 0.375, 0.4, 0.45, 0.475, 0.5]
    pixel3TCAD=[]
    pixel4TCAD=[]
    totalCharge=0

    mg=TMultiGraph("TCAD currents", "TCAD currents")
    leg=TLegend(0.5, 0.5, 0.8, 0.85)
    for p in range(0, len(TCADpositions)):
        for i in range (3, 5):
            filename="/afs/cern.ch/work/n/nalipour/testBeamAnalysis/scripts/TCAD_data/200um_Vbias"+bias+"/eCurrent_strip"+str(i)+"_pos"+str(TCADpositions[p])+".csv"
            gr=returnTCAD(filename, colors[i])
            integralCharge=gr.Integral(0, -1)/echarge
            
            if (i==3):
                pixel3TCAD.append(integralCharge)
            elif (i==4):
                pixel4TCAD.append(integralCharge)

    gr3TCAD=TGraph(len(TCADpositions), array('d', TCADpositions), array('d', pixel3TCAD))
    gr4TCAD=TGraph(len(TCADpositions), array('d', TCADpositions), array('d', pixel4TCAD))

    gr3TCAD.SetMarkerColor(colorPixel3)
    gr4TCAD.SetMarkerColor(colorPixel4)
    
    gr3TCAD.SetMarkerStyle(25)
    gr4TCAD.SetMarkerStyle(25)

    gr3TCAD.SetLineColor(colorPixel3)
    gr4TCAD.SetLineColor(colorPixel4)

    gr3TCAD.SetLineStyle(2)
    gr4TCAD.SetLineStyle(2)

    return gr3TCAD, gr4TCAD

def Theory_chargeInPixel3_4(device, thresh, colorPixel3, colorPixel4):
    pixel3=[]
    pixel4=[]
    posVec=[]
    Threshold=[]
    for i in range (0, 51):
        xpos=i/100.*device.pitch
        chargeVec=device.chargeInEachPixel(80, xpos, "constMob")

        posVec.append(i/100.)
        pixel3.append(chargeVec[2])
        pixel4.append(chargeVec[3])
        Threshold.append(thresh)

    gr3=TGraph(len(posVec), array('d', posVec), array('d', pixel3))
    gr4=TGraph(len(posVec), array('d', posVec), array('d', pixel4))
    grThreshold=TGraph(len(posVec), array('d', posVec), array('d', Threshold))
    grThreshold.SetLineWidth(3)
    grThreshold.SetLineStyle(9)
    grThreshold.SetLineColor(kOrange+2)

    gr3.SetLineColor(colorPixel3)
    gr3.SetMarkerColor(colorPixel3)
    gr3.SetMarkerStyle(31)
    gr4.SetLineColor(colorPixel4)
    gr4.SetMarkerColor(colorPixel4)
    gr4.SetMarkerStyle(31)

    return grThreshold, gr3, gr4

def usage():

    print 'Usage:\n python %s <THL in electrons>' % ( sys.argv[0] )

if __name__ == '__main__':

    gROOT.ProcessLine(".L ~/CLICdpStyle/rootstyle/CLICdpStyle.C");
    CLICdpStyle();

    if (len(sys.argv)<2 or len(sys.argv)>2):
        usage()
        sys.exit( 1 )

    thresh=int(sys.argv[1])

    colors=[kGreen+2, kRed, kBlue, kBlack, kViolet, kMagenta, kOrange, kAzure, kPink]

    # TCAD
    gr3TCAD_35, gr4TCAD_35=TCAD_chargeInPixel3_4(kViolet-1, kRed, "35")
    gr3TCAD, gr4TCAD=TCAD_chargeInPixel3_4(kGreen+2, kBlue, "50")



    # Theory
    device35=Device("n", "p", 200e-4, 35, 30.31, 435, 11.58, 390.6)
    device50=Device("n", "p", 200e-4, 50, 30.31, 435, 11.58, 390.6)
    grThreshold, grTH3_35, grTH4_35=Theory_chargeInPixel3_4(device35, thresh, kViolet-1, kRed)
    grThreshold, grTH3_50, grTH4_50=Theory_chargeInPixel3_4(device50, thresh, kGreen+2, kBlue)

    # Draw 
    mgPixels=TMultiGraph("Charge in different pixels", "Charge in different pixels")
    legPixels=TLegend(0.2, 0.35, 0.5, 0.7)
    mgPixels.Add(gr3TCAD)
    mgPixels.Add(gr4TCAD)
    mgPixels.Add(gr3TCAD_35)
    mgPixels.Add(gr4TCAD_35)
    mgPixels.Add(grTH3_35)
    mgPixels.Add(grTH4_35)
    mgPixels.Add(grTH3_50)
    mgPixels.Add(grTH4_50)
    mgPixels.Add(grThreshold, "l")

    legPixels.AddEntry(gr3TCAD_35, "TCAD: pixel 3: 35 V", "lp")
    legPixels.AddEntry(grTH3_35, "Theory: pixel 3: 35 V", "lp")
    legPixels.AddEntry(gr3TCAD, "TCAD: pixel 3: 50 V", "lp")
    legPixels.AddEntry(grTH3_50, "Theory: pixel 3: 50 V", "lp")
    legPixels.AddEntry(gr4TCAD_35, "TCAD: pixel 4: 35 V", "lp")
    legPixels.AddEntry(grTH4_35, "Theory: pixel 4: 35 V", "lp")
    legPixels.AddEntry(gr4TCAD, "TCAD: pixel 4: 50 V", "lp")
    legPixels.AddEntry(grTH4_50, "Theory: pixel 4: 50 V", "lp")

    mgPixels.Draw("ALP")
    legPixels.Draw()
    mgPixels.GetXaxis().SetTitle("Hit position X/pitch")
    mgPixels.GetYaxis().SetTitle("Collected charge [e-]")
    mgPixels.GetYaxis().SetLabelSize(0.05)

    gPad.Update()
    gPad.Modified()

    bla=raw_input()
    
