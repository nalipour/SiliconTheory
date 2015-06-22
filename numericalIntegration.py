import sys

from ROOT import *
import ROOT
from array import array
from constants_theory import *
from theory import *


def ReturnDriftTime(device):

    zpos=[]

    ConstMob=[]
    NonConstMob=[]

    func_ConstMob=[]
    func_NonConstMob=[]
    
    driftTimeTheory=[]
    E=[]

    for z in range (0, int(device.thickness*1e4)+1):
        z=z/10000.0 #cm
        zpos.append(z)
            
        Ef, mob, t_d=device.Efield(z)
        driftTimeTheory.append(t_d)
        
        ConstMob.append(device.CarrierMobilityConst)
        NonConstMob.append(mob)
        
        func_ConstMob.append(1./(device.CarrierMobilityConst*Ef))
        func_NonConstMob.append(1./(mob*Ef))
        E.append(Ef)
        
    grDriftTimeTheory=TGraph(len(zpos), array('f', zpos), array('f', driftTimeTheory))
    grE=TGraph(len(zpos), array('f', zpos), array('f', E))

    grConstMob=TGraph(len(zpos), array('f', zpos), array('f', ConstMob))
    grNonConstMob=TGraph(len(zpos), array('f', zpos), array('f', NonConstMob))
    grConstMob.SetLineColor(kGreen+2)

    grFuncConstMob=TGraph(len(zpos), array('f', zpos), array('f', func_ConstMob)) 
    grFuncNonConstMob=TGraph(len(zpos), array('f', zpos), array('f', func_NonConstMob)) 
    grFuncConstMob.SetLineColor(kGreen+2)
    

    # print "Integral constMob=", TMath.Sqrt(2*device.DiffusionConst*grFuncConstMob.Integral(0, -1))*1e4
    print "Integral constMob=", NumInt(grFuncConstMob, 0, 0.02, 10000) #NOT CORRECT -> WHYYYYYY????
    print "Integral NonConstMob=", NumInt(grFuncNonConstMob, 0, 0.02, 10000)
    print "Integral const mu=", NumInt(grConstMob, 0, 0.02, 10000)

    return grE, grConstMob, grNonConstMob, grConstMob, grFuncConstMob, grFuncNonConstMob, grDriftTimeTheory


def IntegrateGraph(gr, device, Nstep):

    IntegralTotal=[]
    zpos=[]
    for z in range (0, int(device.thickness*1e4)+1):
        z=z/10000.0 #cm
        zpos.append(z)
        temp=NumInt(gr, 0.0, z, Nstep)

        IntegralTotal.append(temp)

#    print "size IntegralTotal=", len(IntegralTotal)
    grIntegral=TGraph(len(zpos), array('f', zpos), array('f', IntegralTotal))
    return grIntegral


def NumInt(gr, x_min, x_max, Nstep):
    gr_min=TMath.MinElement(gr.GetN(), gr.GetX())
    gr_max=TMath.MaxElement(gr.GetN(), gr.GetX())

    
    if (x_min<gr_min):
        x_min=gr_min
    if(x_max>gr_max):
        x_max=gr_max

    stepSize=(x_max-x_min)/Nstep
    Integral=0
    for i in range(0, int(Nstep)):
        x1=x_min+i*stepSize
        x2=x_min+(i+1)*stepSize

        y1=gr.Eval(x1)
        y2=gr.Eval(x2)
        Integral+=(y1+y2)*(x2-x1)/2.

    return Integral



if __name__ == '__main__':

    gROOT.ProcessLine(".L ~/CLICdpStyle/rootstyle/CLICdpStyle.C")
    CLICdpStyle()
    
    
    device=Device("n", "p", 200e-4, 50, 30.31, 435, 11.58, 390.6) #B06

    grE, grConstMob, grNonConstMob, grConstMob, grFuncConstMob, grFuncNonConstMob, grDriftTimeTheory=ReturnDriftTime(device)


    # Integrate time
    grTimeConstMob=IntegrateGraph(grFuncConstMob, device, 10000)
    grTimeNonConstMob=IntegrateGraph(grFuncNonConstMob, device, 10000)

    grDriftTimeTheory.SetLineColor(kRed)
    grDriftTimeTheory.SetMarkerStyle(29)
    grDriftTimeTheory.SetMarkerColor(kRed)
    grTimeConstMob.SetLineColor(kGreen+2)
    grTimeNonConstMob.SetLineColor(kBlue)

    mg=TMultiGraph("mg", "mg")
    leg=TLegend(0.2, 0.5, 0.5, 0.7)

    mg.Add(grDriftTimeTheory)
    leg.AddEntry(grDriftTimeTheory, "Theory", "p")
    mg.Add(grTimeConstMob)
    leg.AddEntry(grTimeConstMob, "Const mob", "l")
    mg.Add(grTimeNonConstMob)
    leg.AddEntry(grTimeNonConstMob, "Non const mob", "l")

    canv=TCanvas("can", "can")
    mg.Draw("ALP")
    leg.Draw()
    mg.GetXaxis().SetTitle("Depth [cm]")
    mg.GetYaxis().SetTitle("Drift time [s]")
    mg.GetXaxis().SetRangeUser(0, device.thickness*1e4)
    canv.Update()
    canv.Print("plots/driftTime.pdf")

    # Draw
    canvE=TCanvas("Efield", "Efield")
    grE.Draw("ALP")
    grE.GetXaxis().SetTitle("Depth [cm]")
    grE.GetYaxis().SetTitle("|E_{field}| [V/cm]")
    grE.GetXaxis().SetRangeUser(0, device.thickness*1e4)
    canvE.Update()
    canvE.Print("plots/Efield.pdf")

    canvMob=TCanvas("Mobility", "Mobility")
    mgMob=TMultiGraph("mob", "mob")
    mgMob.Add(grConstMob)
    mgMob.Add(grNonConstMob)
    legMob=TLegend(0.55, 0.2, 0.8, 0.4)
    legMob.AddEntry(grConstMob, "#mu_{const}", "l")
    legMob.AddEntry(grNonConstMob, "#mu_{non-constant}", "l")

    mgMob.Draw("ALP")
    legMob.Draw()
    mgMob.GetXaxis().SetTitle("Depth [cm]")
    mgMob.GetYaxis().SetTitle("#mu_{carrier} [cm^{2}/Vs]")
    mgMob.GetXaxis().SetRangeUser(0, device.thickness*1e4)
    canvMob.Update()
    canvMob.Print("plots/Mobility.pdf")


    canvFuncMob=TCanvas("Func Mob", "Func Mob")
    mgFunc=TMultiGraph("Func mob", "Func mob")
    mgFunc.Add(grFuncConstMob)
    mgFunc.Add(grFuncNonConstMob)
    legFunc=TLegend(0.2, 0.5, 0.5, 0.7)
    legFunc.AddEntry(grFuncConstMob, "#mu_{const}", "l")
    legFunc.AddEntry(grFuncNonConstMob, "#mu_{non-constant}", "l")

    mgFunc.Draw("ALP")
    legFunc.Draw()
    mgFunc.GetXaxis().SetTitle("Depth [cm]")
    mgFunc.GetYaxis().SetTitle("1/(#mu(z)E(z)) [s/cm]")
    mgFunc.GetXaxis().SetRangeUser(0, device.thickness*1e4)
    canvFuncMob.Update()
    canvFuncMob.Print("plots/FunctionToIntegrate.pdf")


    # Calculation of the diffusion
    zIntegral=[]
    timeIntegrated=[]
    timeLin=[]
    diffusionConst=[]
    diffusionConstTheory=[]
    sigmaVarMob=[]
    sigma=[]
    for depth in range(0, int(device.thickness*10000.0)+1):
        zpos=[]
        E=[]
        mobility=[]
        t_drift=[]
        func=[]
        zIntegral.append(depth)
        Ef, mob, t_d=device.Efield(depth/10000.0)

        timeLin.append(t_d)

        diffusion_val=TMath.K()*300*mob/echarge
        timeIntegrated_val=0
        diffusionConst.append(diffusion_val)
        diffusionConstTheory.append(device.DiffusionConst)
        

#        print grTimeNonConstMob.Eval(depth/10000.)
        sigmaVarMob.append(TMath.Sqrt(2*diffusion_val*grTimeNonConstMob.Eval(depth/10000.))*10000)
        sigma.append(TMath.Sqrt(2*device.DiffusionConst*grTimeConstMob.Eval(depth/10000.))*10000)

        


    
    # ===== Draw Diffusion Constant ===== #
    canvDiffusion=TCanvas("Diffusion", "Diffusion")
    grDiffusionConst=TGraph(len(zIntegral), array('f', zIntegral), array('f', diffusionConst))
    grDiffusionConstTheory=TGraph(len(zIntegral), array('f', zIntegral), array('f', diffusionConstTheory))

    grDiffusionConst.SetMarkerColor(kBlue)
    grDiffusionConst.SetMarkerStyle(7)
    grDiffusionConst.SetLineColor(kBlue)

    grDiffusionConstTheory.SetMarkerColor(kBlack)
    grDiffusionConstTheory.SetMarkerStyle(7)
    grDiffusionConstTheory.SetLineColor(kBlack)

    mgDiffConst=TMultiGraph("diffusion constant", "diffusion constant")
    mgDiffConst.Add(grDiffusionConst, "l")
    mgDiffConst.Add(grDiffusionConstTheory, "l")
    legDiff=TLegend(0.2, 0.5, 0.5, 0.7)
    legDiff.AddEntry(grDiffusionConst, "#mu_{non-constant}", "l")
    legDiff.AddEntry(grDiffusionConstTheory,  "#mu_{const}", "l")
    mgDiffConst.Draw("ALP")
    legDiff.Draw()
    mgDiffConst.GetXaxis().SetTitle("Depth [#mum]")
    mgDiffConst.GetYaxis().SetTitle("Diffusion Constant [cm^{2}s^{-1}]")
    mgDiffConst.GetXaxis().SetRangeUser(0, device.thickness*1e4)
    canvDiffusion.Update()
    canvDiffusion.Modified()
    canvDiffusion.Print("plots/DiffusionConstant.pdf")
    # ===== End Draw Diffusion Constant ===== #

    # ===== Draw Diffusion sigma ===== #
    canvSigma=TCanvas("sigma", "sigma")
    grSigmaVar=TGraph(len(zIntegral), array('f', zIntegral), array('f', sigmaVarMob))
    grSigma=TGraph(len(zIntegral), array('f', zIntegral), array('f', sigma))
    
    grSigma.SetMarkerColor(kBlack)
    grSigma.SetMarkerStyle(7)
    grSigma.SetLineColor(kBlack)

    grSigmaVar.SetMarkerColor(kBlue)
    grSigmaVar.SetMarkerStyle(7)
    grSigmaVar.SetLineColor(kBlue)

    mgSigma=TMultiGraph("Sigma", "Sigma")
    mgSigma.Add(grSigma, "l")
    mgSigma.Add(grSigmaVar, "l")
    legSigma=TLegend(0.2, 0.6, 0.5, 0.8)
    legSigma.AddEntry(grSigmaVar, "#mu_{non-constant}", "l")
    legSigma.AddEntry(grSigma,  "#mu_{const}", "l")
    mgSigma.Draw("ALP")
    legSigma.Draw()
    mgSigma.GetXaxis().SetTitle("Depth [#mum]")
    mgSigma.GetYaxis().SetTitle("#sigma_{diffusion} [#mum]")
    mgSigma.GetXaxis().SetRangeUser(0, device.thickness*1e4)
    canvSigma.Update()
    canvSigma.Modified()
    canvSigma.Print("plots/SigmaDiffusion.pdf")
    # ===== End Draw Diffusion sigma ===== #


    bla=raw_input()
