import sys

from ROOT import *
import ROOT
from array import array
from constants_theory import *
from theory import *


def calculations(biasOrTHLScan, thresholds, biasVoltages, carrierType, bulkType, thickness, V_D, thresholdTrans_a, thresholdTrans_b):
    colorsVec=[kGreen+2, kRed, kBlue, kBlack, kViolet, kMagenta, kOrange, kAzure, kPink, kYellow, kGray, kTeal, kCyan, kViolet-1, kAzure+10, kCyan+4]
    ThresholdInElectrons=[]
    ThresholdCrossing=[]
    mgAll=TMultiGraph("mg", "mg")
    for b in range(0, len(biasVoltages)):
        for t in range(0, len(thresholds)):
            device=Device(carrierType, bulkType, thickness, biasVoltages[b], V_D, thresholds[t], thresholdTrans_a, thresholdTrans_b)
            pixel=[]
            posVec=[]
            for i in range (0, 51):
                xpos=i/100.*device.pitch
                chargeVec=device.chargeInEachPixel(80, xpos, "constMob")
                posVec.append(i/100.)
                pixel.append(chargeVec[3])
                
            # device.threshold_electrons=device.threshold_electrons-800
            gr=TGraph(len(posVec), array('d', posVec), array('d', pixel))
            gr.SetLineColor(colorsVec[b*len(thresholds)+t])
            fit=TF1("FitPixel", "pol15", 0, 0.5)
            fit.SetLineColor(colorsVec[b*len(thresholds)+t])
            gr.Fit("FitPixel", "R")
            print "Threshold [e] = ", device.threshold_electrons
            ThresholdInElectrons.append(device.threshold_electrons)

            mgAll.Add(gr)

            xThresh=fit.GetX(device.threshold_electrons)
            ThresholdCrossing.append(xThresh)

            ratioSize1=[]

            for i in range(0, len(ThresholdCrossing)):
                fraction_1Hit=((ThresholdCrossing[i]*2)**2)
                ratioSize1.append(fraction_1Hit)
                
    canv=TCanvas("Integral", "Integral")
    if (biasOrTHLScan=="bias"):
        
        grSize1=TGraph(len(biasVoltages), array('d', biasVoltages), array('d', ratioSize1))
        grSize1.SetMarkerStyle(24)
        grSize1.SetMarkerColor(colors)
        grSize1.SetLineColor(colors)
        return grSize1, mgAll
    
    elif (biasOrTHLScan=="THL"):
        # grSize1=TGraph(len(thresholds), array('d', thresholds), array('d', ratioSize1))
        grSize1=TGraph(len(ThresholdInElectrons), array('d', ThresholdInElectrons), array('d', ratioSize1))
        grSize1.SetMarkerStyle(24)
        grSize1.SetMarkerColor(colors)
        grSize1.SetLineColor(colors)
        return grSize1, mgAll, ThresholdInElectrons


def clusterSize(runs, xlabels, color, changeLabels):

    ratios1hit=[]
    ratios2hit=[]
    index=[]
    for i in range(0, len(runs)):
        index.append(i)
        runstr=str(runs[i])
        filename="/afs/cern.ch/work/s/sredford/pixel_pyeudetdata/Run"+runstr+"/EtaCorrection/pyEudetNtuple_run"+runstr+"_EtaCorrection.root"
        file=TFile.Open(filename)
        tree=file.Get("clusters")
        
        canv=TCanvas("clusize", "clusize")
        hist=TH1D("hist", "hist", 10, 0, 10)
        tree.Draw("size >> hist")
        hist.Scale(1/hist.Integral())
        
        hist.GetXaxis().SetTitle("Cluster size");
        hist.GetYaxis().SetTitle("Entries");
        # canv.Print("plots/clusterSize/Run"+runstr+".pdf")
        
        ratio1=hist.GetBinContent(2)    # 1-hit clusters
        ratio2=hist.GetBinContent(3)    # 2-hit clusters

        ratios1hit.append(ratio1)
        ratios2hit.append(ratio2)


    if (changeLabels!=""):
        gr1hit=TGraph(len(index), array('d', index), array('d', ratios1hit))
        gr1hit.SetMarkerStyle(25)
        gr1hit.SetMarkerColor(color)
        
        
        xaxis=gr1hit.GetXaxis()
        for i in range(0, len(ratios1hit)):
            binIndex = gr1hit.GetXaxis().FindBin(i)
            gr1hit.GetXaxis().SetBinLabel(binIndex, str(xlabels[i]))

 
        return ratios1hit, ratios2hit, gr1hit
    else:
        gr1hit=TGraph(len(xlabels), array('d', xlabels), array('d', ratios1hit))
        gr1hit.SetMarkerStyle(25)
        gr1hit.SetMarkerColor(color)
        gr1hit.SetLineColor(color)

        gr2hit=TGraph(len(xlabels), array('d', xlabels), array('d', ratios2hit))
        gr2hit.SetMarkerStyle(24)
        gr2hit.SetMarkerColor(kRed)
        gr2hit.SetLineColor(kRed)
        return ratios1hit, ratios2hit, gr1hit, gr2hit


if __name__ == '__main__':

    gROOT.ProcessLine(".L ~/CLICdpStyle/rootstyle/CLICdpStyle.C");
    CLICdpStyle();

    # ==== B06 ==== #
    # Normal runs
    runs=[2310, 2311, 2311, 2312, 2357, 2358, 2359, 2360, 2361, 2362, 2363, 2364] # Normal runs 
    ratios1hit, ratios2hit, gr1hit=clusterSize(runs, runs, kBlack, "No")
    
    # Bias scan
    runs_bias=[2296, 2297, 2298, 2300, 2301, 2302, 2303, 2304, 2305, 2306, 2307, 2308, 2309]
    runs_voltages=[5, 10, 15, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]
    # ratios1hit_bias, ratios2hit_bias, gr1hit_bias, gr2hit_bias=clusterSize(runs_bias, runs_voltages, kBlack, "")

    # THL scan
    runs_threshold=[2317, 2319, 2320, 2321, 2322, 2323, 2324, 2325, 2326, 2327, 2328, 2329, 2330, 2331, 2333, 2334]
    runs_THL=[425, 435, 440, 445, 450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 505, 510]
    # ratios1hit_THL, ratios2hit_THL, gr1hit_THL, gr2hit_THL=clusterSize(runs_threshold, runs_THL, kBlack, "")

    # Parameters for theory
    V_D=30.31
    V_B=50
    threshold=1066
    thickness=200e-4
    temperature=300
    charge=80
    colors=kBlue
    z=200e-4
    thresholdTrans_a=11.58
    thresholdTrans_b=390.6
    opTHL=435
    carrierType="n"
    bulkType="p"


    # # ==== A06 ==== # -------> No normal runs in  pyEudet
    # # # Normal runs
    # # runs=[153] # Normal runs 
    # # ratios1hit, ratios2hit, gr1hit=clusterSize(runs, runs, kBlack, "No")
    
    # # Bias scan
    # runs_bias=[127, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141]
    # runs_voltages=[15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45]
    # ratios1hit_bias, ratios2hit_bias, gr1hit_bias, gr2hit_bias=clusterSize(runs_bias, runs_voltages, kBlack, "")

    # # THL scan
    # runs_threshold=[1048, 1049, 1053, 1054, 1057, 1058, 1060, 1061, 1062]
    # runs_THL=[322, 324, 324, 326, 328, 330, 334, 336, 338]
    # ratios1hit_THL, ratios2hit_THL, gr1hit_THL, gr2hit_THL=clusterSize(runs_threshold, runs_THL, kBlack, "")

    # # Parameters for theory
    # V_D=10
    # V_B=15
    # threshold=855
    # thickness=50e-4
    # temperature=300
    # charge=80
    # colors=kBlue
    # z=50e-4
    # thresholdTrans_a=-12.36
    # thresholdTrans_b=364
    # opTHL=326
    # carrierType="p"
    # bulkType="n"

#     # ==== L04 ==== #
#     # Bias scan
#     runs_bias=[220, 221, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235]
#     runs_voltages=[5, 7.5, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5]
#     # THL scan
#     runs_threshold=[1174, 1175, 1177, 1178, 1190, 1180, 1189]#, 1191, 1192, 1193]
#     runs_THL=[385, 390, 400, 405, 410, 415, 435]#, 410, 410, 410, 410]
#     # Parameters for theory
#     V_D=19.64
#     V_B=35
# #    threshold=855
#     thickness=100e-4
#     temperature=300
#     charge=80
#     colors=kBlue
#     z=100e-4
#     thresholdTrans_a=-11.68
#     thresholdTrans_b=448.5
#     opTHL=410
#     carrierType="p"
#     bulkType="n"
#     # ======= Data ======= #
    
    ratios1hit_bias, ratios2hit_bias, gr1hit_bias, gr2hit_bias=clusterSize(runs_bias, runs_voltages, kBlack, "")
    ratios1hit_THL, ratios2hit_THL, gr1hit_THL, gr2hit_THL=clusterSize(runs_threshold, runs_THL, kBlack, "")
    # bla=raw_input()
    # ===== Theory ===== #
    # Theory bias scan
    # theory_biasScan_size1, V_Bs, ThresholdCrossing=compareDifferentBiasVoltage(V_D, runs_voltages, thickness, temperature, 80, colors, threshold)

    # Theory THL scan
    # theory_THLscan_size1=theorycalulations_THLscan(V_D, V_B, thickness, temperature, charge, colors, runs_THL, a, b)
    theory_THLscan_size1, mgAllTHL, ThresholdInElectrons=calculations("THL", runs_THL, [V_B], carrierType, bulkType, thickness, V_D, thresholdTrans_a, thresholdTrans_b)
    theory_biasScan_size1, mgAllBias=calculations("bias", [opTHL], runs_voltages, carrierType, bulkType, thickness, V_D, thresholdTrans_a, thresholdTrans_b)

    
    
    # # Draw normal
    # canvNormal=TCanvas("normal runs", "normal runs")
    # gr1hit.Draw("AP")
    # gPad.Modified()
    # gPad.Update()


    # Draw Bias scan
    canvNormalBias=TCanvas("bias runs", "bias runs")
    legBias=TLegend(0.55, 0.2, 0.8, 0.4)
    mgBias=TMultiGraph("Bias", "Bias")
    mgBias.Add(gr1hit_bias)
    mgBias.Add(theory_biasScan_size1)
    
    legBias.AddEntry(gr1hit_bias, "Data", "lp")
    legBias.AddEntry(theory_biasScan_size1, "Theory", "lp")
    

    mgBias.Draw("ALP")
    legBias.Draw()
    legBias.SetFillStyle(0)
    mgBias.GetXaxis().SetTitle("|V_{bias}| [V]")
    mgBias.GetYaxis().SetTitle("Fraction of 1-hit clusters")
    gPad.Modified()
    gPad.Update()

    #Draw THL scan
    dataTHL1=TGraph(len(ThresholdInElectrons), array('d', ThresholdInElectrons), array('d', ratios1hit_THL))
    dataTHL1.SetMarkerStyle(25)
    dataTHL1.SetMarkerColor(kBlack)
    dataTHL1.SetLineColor(kBlack)

    print "ThresholdInElectrons=", ThresholdInElectrons
    canvNormalTHL=TCanvas("THL runs", "THL runs")
    legTHL=TLegend(0.55, 0.2, 0.8, 0.4)
    mgTHL=TMultiGraph("THL", "THL")
    mgTHL.Add(dataTHL1)
    # mgTHL.Add(gr1hit_THL)
    # mgTHL.Add(gr2hit_THL)
    mgTHL.Add(theory_THLscan_size1)

    legTHL.AddEntry(dataTHL1, "Data",  "lp")
    # legTHL.AddEntry(gr1hit_THL, "1-hit",  "lp")
#    legTHL.AddEntry(gr2hit_THL, "2-hit", "lp")
    legTHL.AddEntry(theory_THLscan_size1, "Theory", "lp")

    mgTHL.Draw("ALP")
    legTHL.Draw()
    legTHL.SetFillStyle(0)
    mgTHL.GetXaxis().SetTitle("Threshold [e-]")
    mgTHL.GetYaxis().SetTitle("Fraction of 1-hit clusters")
    gPad.Modified()
    gPad.Update()

    bla=raw_input()
