import sys

from ROOT import *
import ROOT
from array import array
from constants_theory import *
from theory import *


def FunctionToIntegrate(device):

    zpos=[]

    ConstMob=[]
    NonConstMob=[]

    func_ConstMob=[]
    func_NonConstMob=[]
    
    E=[]

    for z in range (0, int(device.thickness*1e4)+1):
        z=z/10000.0 #cm
        zpos.append(z)
            
        Ef, mob, t_d=device.Efield(z)
        # print z, "=======", t_d
        
        ConstMob.append(device.CarrierMobilityConst)
        NonConstMob.append(mob)
        
        func_ConstMob.append(1./(device.CarrierMobilityConst*Ef))
        func_NonConstMob.append(1./(mob*Ef))
        E.append(Ef)
        
    grE=TGraph(len(zpos), array('f', zpos), array('f', E))

    grConstMob=TGraph(len(zpos), array('f', zpos), array('f', ConstMob))
    grNonConstMob=TGraph(len(zpos), array('f', zpos), array('f', NonConstMob))
    grConstMob.SetLineColor(kGreen+2)

    grFuncConstMob=TGraph(len(zpos), array('f', zpos), array('f', func_ConstMob)) 
    grFuncNonConstMob=TGraph(len(zpos), array('f', zpos), array('f', func_NonConstMob)) 
    grFuncConstMob.SetLineColor(kGreen+2)
    
    # print "Integral constMob=", TMath.Sqrt(2*device.DiffusionConst*grFuncConstMob.Integral(0, -1))*1e4
    # print "Integral constMob=", NumInt(grFuncConstMob, 0, 0.02, 10000) 
    # print "Integral NonConstMob=", NumInt(grFuncNonConstMob, 0, 0.02, 10000)
    # print "Integral const mu=", NumInt(grConstMob, 0, 0.02, 10000)

    # Draw
    canvE=TCanvas("Efield", "Efield")
    grE.Draw("ALP")
    grE.GetXaxis().SetTitle("Depth [#mum]")
    grE.GetYaxis().SetTitle("|E_{field}| [V/cm]")
    grE.GetXaxis().SetRangeUser(0, device.thickness*1e4)
    canvE.Update()

    canvMob=TCanvas("Mobility", "Mobility")
    mgMob=TMultiGraph("mob", "mob")
    mgMob.Add(grConstMob)
    mgMob.Add(grNonConstMob)
    legMob=TLegend(0.55, 0.2, 0.8, 0.4)
    legMob.AddEntry(grConstMob, "#mu_{const}", "l")
    legMob.AddEntry(grNonConstMob, "#mu_{non-constant}", "l")

    mgMob.Draw("ALP")
    legMob.Draw()
    mgMob.GetXaxis().SetTitle("Depth [#mum]")
    mgMob.GetYaxis().SetTitle("#mu_{carrier} [cm^{2}/Vs]")
    mgMob.GetXaxis().SetRangeUser(0, device.thickness*1e4)
    canvMob.Update()


    canvFuncMob=TCanvas("Func Mob", "Func Mob")
    mgFunc=TMultiGraph("Func mob", "Func mob")
    mgFunc.Add(grFuncConstMob)
    mgFunc.Add(grFuncNonConstMob)
    legFunc=TLegend(0.55, 0.2, 0.8, 0.4)
    legFunc.AddEntry(grFuncConstMob, "#mu_{const}", "l")
    legFunc.AddEntry(grFuncNonConstMob, "#mu_{non-constant}", "l")

    mgFunc.Draw("ALP")
    legFunc.Draw()
    mgFunc.GetXaxis().SetTitle("Depth [#mum]")
    mgFunc.GetYaxis().SetTitle("1/(#muE) [s/cm]")
    mgFunc.GetXaxis().SetRangeUser(0, device.thickness*1e4)
    canvFuncMob.Update()
    

    return 


def NumInt(gr, x_min, x_max, Nstep):
 
    gr_min=TMath.MinElement(gr.GetN(), gr.GetX())
    gr_max=TMath.MaxElement(gr.GetN(), gr.GetX())
    
    if (x_min<gr_min):
        x_min=gr_min
    if(x_max>gr_max):
        x_max=gr_max

    stepSize=(x_max-x_min)/Nstep

    Integral=0
    for i in range(0, Nstep):
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

    FunctionToIntegrate(device)


    zIntegral=[]
    timeIntegrated=[]
    timeLin=[]
    diffusionConst=[]
    sigmaVarMob=[]
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
        
        
        for z in range (0, depth):
            z=z/10000.0 #cm
            zpos.append(z)
            
            Ef, mob, t_d=device.Efield(z)
            func.append(1./(mob*Ef))
            E.append(Ef)
            mobility.append(mob)
    
        if (zpos!=[]):
            gr=TGraph(len(zpos), array('f', zpos), array('f', func))
            gr.SetMarkerColor(kRed)
            gr.SetMarkerStyle(7)
            gr.SetLineColor(kRed)

            # gr.Draw("ALP")
            timeIntegrated_val=gr.Integral(0, -1)
            

        timeIntegrated.append(timeIntegrated_val)
        sigmaVarMob.append(TMath.Sqrt(2*diffusion_val*timeIntegrated_val)*10000)
    
    grLin=TGraph(len(zIntegral), array('f', zIntegral), array('f', timeLin))
    grIntegrated=TGraph(len(zIntegral), array('f', zIntegral), array('f', timeIntegrated))



    

    grLin.SetMarkerColor(kRed)
    grLin.SetMarkerStyle(7)
    grLin.SetLineColor(kRed)

    grIntegrated.SetMarkerColor(kBlue)
    grIntegrated.SetMarkerStyle(7)
    grIntegrated.SetLineColor(kBlue)

    canvTime=TCanvas("Time", "Time")
 
    grIntegrated.Draw("ALP")
    grLin.Draw("same")



    canvDiffusion=TCanvas("Diffusion", "Diffusion")
    grDiffusionConst=TGraph(len(zIntegral), array('f', zIntegral), array('f', diffusionConst))
    grDiffusionConst.SetMarkerColor(kBlue)
    grDiffusionConst.SetMarkerStyle(7)
    grDiffusionConst.SetLineColor(kBlue)
    grDiffusionConst.Draw("ALP")


    canvSigma=TCanvas("sigma", "sigma")
    grSigma=TGraph(len(zIntegral), array('f', zIntegral), array('f', sigmaVarMob))
    grSigma.SetMarkerColor(kBlue)
    grSigma.SetMarkerStyle(7)
    grSigma.SetLineColor(kBlue)

    grE, grM, grS, grS_varMob=device.returnArray(kBlue, kGreen+2, kBlack)

    grSigma.Draw("ALP")
    grS.Draw("same")


    bla=raw_input()
