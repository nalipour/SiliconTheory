import sys

from ROOT import *
import ROOT
from array import array
from constants_theory import *


class Device:
    "Device properties"
    
    carrierType=""
    bulkType=""

    thickness=0 # [cm]
    V_B=0
    V_D=0
    pitch=55 # [um]

    temperature=300
    thresholdDAC=0
    thresholdTrans_a=0
    thresholdTrans_b=0
    

    vm=0
    Ec=0
    beta=0
    DiffusionConst=0
    CarrierMobilityConst=0


    
    def __init__(self, carrierType, bulkType, thickness, V_B, V_D, thresholdDAC, thresholdTrans_a, thresholdTrans_b):

        self.carrierType=carrierType
        self.bulkType=bulkType
        self.thickness=thickness
        self.V_B=V_B
        self.V_D=V_D

        self.thresholdTrans_a=thresholdTrans_a
        self.thresholdTrans_b=thresholdTrans_b
        self.threshold_electrons=((thresholdDAC-self.thresholdTrans_b)/self.thresholdTrans_a)/3.6*1000.0
        
        if (self.carrierType=="p"):
            self.vm=1.62*TMath.Power(10, 8)*TMath.Power(self.temperature, -0.52)
            self.Ec=1.24*TMath.Power(self.temperature, 1.68)
            self.beta=0.46*TMath.Power(self.temperature, 0.17)
            self.DiffusionConst=Default_Hole_D
            self.CarrierMobilityConst=Default_Hole_Mobility

        elif (self.carrierType=="n"):
            self.vm=1.53*TMath.Power(10, 9)*TMath.Power(self.temperature, -0.87)
            self.Ec=1.01*TMath.Power(self.temperature, 1.55)
            self.beta=2.57*TMath.Power(10, -2)*TMath.Power(self.temperature, 0.66)
            self.DiffusionConst=Default_Electron_D
            self.CarrierMobilityConst=Default_Electron_Mobility


    def Efield(self, z):
        Efield=((self.V_B-self.V_D)/self.thickness+(1-z/self.thickness)*2*self.V_D/self.thickness) #[V/cm]
        mobility_c=(self.vm/self.Ec)/(TMath.Power((1+TMath.Power(TMath.Abs(Efield)/self.Ec, self.beta)), 1.0/self.beta))
        # t_drift=(self.thickness*self.thickness)/(2*mobility_c*self.V_D)*TMath.Log((self.V_B+self.V_D)/(self.V_B+self.V_D-2.0*self.V_D*z/self.thickness))
        t_drift=(self.thickness*self.thickness)/(2*self.CarrierMobilityConst*self.V_D)*TMath.Log((self.V_B+self.V_D)/(self.V_B+self.V_D-2.0*self.V_D*z/self.thickness))
        return Efield, mobility_c, t_drift

    def sigmaDiffusion(self, z):
        # Calculates the sigma of diffusion
        # V_B: bias voltage [V]
        # V_D: depletion voltage [V]
        # thickness of the sensor [cm]
        # temperature: [K]
        # TMath.K(): Boltzmann constant
        # The diffusion sigma is in [mum]

        # z: Depth [cm] at which sigma is calculated
        depletionDepth=(self.thickness/(2*self.V_D))*(self.V_D+self.V_B)-1e-4
        if (z<depletionDepth):
            sigma=(TMath.Sqrt((TMath.K()*self.temperature*self.thickness*self.thickness/(echarge*self.V_D))*TMath.Log((self.V_B+self.V_D)/TMath.Abs(self.V_B+self.V_D-2.0*self.V_D*z/self.thickness))))*10000. #mum
            return sigma
            #return TMath.Sqrt(4)*sigma
        else:
            return 0
    
    
    def returnArray(self, colorE, colorM, colorS):
        zpos=[]
        E=[]
        sigma=[]
        mobility=[]
        t_drift=[]
        sigma_varMob=[]

        for z in range (0, int(self.thickness*10000.0)+1):
            zpos.append(z)
            z=z/10000.0 #cm
            ef, mob, t_d=self.Efield(z)
            E.append(ef)
            mobility.append(mob)
            t_drift.append(t_d)
            sigma.append(self.sigmaDiffusion(z)) 
            sigma_varMob.append(TMath.Sqrt(2*self.DiffusionConst*t_d*1e8))
            
        gr_Efield=TGraph(len(zpos), array('f', zpos), array('f', E))
        gr_mobility=TGraph(len(zpos), array('f', zpos), array('f', mobility))
        gr_Sigma=TGraph(len(zpos), array('f', zpos), array('f', sigma))
        gr_Sigma_varMob=TGraph(len(zpos), array('f', zpos), array('f', sigma_varMob))        


        gr_Efield.SetMarkerColor(colorE)
        gr_Efield.SetMarkerStyle(7)
        gr_Efield.SetLineColor(colorE)
        
        gr_mobility.SetMarkerColor(colorM)
        gr_mobility.SetMarkerStyle(7)
        gr_mobility.SetLineColor(colorM)
        
        gr_Sigma.SetLineColor(colorS)
        gr_Sigma_varMob.SetLineColor(colorE)

        return gr_Efield, gr_mobility, gr_Sigma, gr_Sigma_varMob

    def integrateGaussian(self, xhit, Sigma, x1, x2):
        Integral=(-TMath.Erf((x1-xhit)/(TMath.Sqrt(2.)*Sigma))+TMath.Erf((x2-xhit)/(TMath.Sqrt(2.)*Sigma)))/2.0
        return Integral

    def chargeInEachPixel(self, charge, xhit, diffusionMethod):
        pixels=[-2, -1, 0, 1, 2]
        chargePerPixel=[0, 0, 0, 0, 0]
        for z in range(0, int(self.thickness*1e4+1)):
            z=z*1e-4
            sigma_z=0
            if (diffusionMethod=="constMob"):
                sigma_z=self.sigmaDiffusion(z)


            if (sigma_z>0):
                for i in range(0, len(pixels)):
                    integral=self.integrateGaussian(xhit, sigma_z, (-self.pitch/2.0 + pixels[i]*self.pitch), (-self.pitch/2.0+(pixels[i]+1)*self.pitch))*charge
                    chargePerPixel[i]+=integral
                    
        return chargePerPixel



if __name__ == '__main__':
    
    gROOT.ProcessLine(".L ~/CLICdpStyle/rootstyle/CLICdpStyle.C");
    CLICdpStyle();
    gStyle.SetPadTickY(0)
    gStyle.SetTitleOffset(1,"x")
    gStyle.SetLabelSize(0.05,"xyz")
    gStyle.SetTitleOffset(1.1,"yz")
    gStyle.SetTextSize(0.03)
    gStyle.SetPadLeftMargin(0.17)
    gStyle.SetPadRightMargin(0.17)
    
#    gStyle.SetMarkerSize(0.06)
# (self, carrierType, bulkType, thickness, V_B, V_D, thresholdDAC, thresholdTrans_a, thresholdTrans_b)
    # device=Device("n", "p", 200e-4, 35, 30.31, 435, 11.58, 390.6) #B06
    # device=Device("p", "n", 50e-4, 15, 10, 326, -12.36, 364)#A06
    device=Device("p", "n", 100e-4, 35, 19.64, 410, -11.68, 448.5)#L04
    


    device.Efield(100e-4)
    grE, grM, grS, grS_varMob=device.returnArray(kBlue, kGreen+2, kBlack)

    canvM=TCanvas("mob", "mob")
    grM.Draw("ALP")


    canvE=TCanvas("Efield", "Efield")
 
    grE.Draw("ALP")
    grE.GetXaxis().SetTitle("Depth [#mum]")
    grE.GetYaxis().SetTitle("|E_{field}| [V/cm]")
    grE.GetXaxis().SetRangeUser(0, device.thickness*1e4)
    canvE.Update()

    rightmax = 1.1*TMath.MaxElement(grM.GetN(), grM.GetY())#grM.GetMaximum()
    rightmin = TMath.MinElement(grM.GetN(), grM.GetY())#grM.GetMaximum()
    print rightmax-rightmin
    print "rightmax=", rightmax
    print "rightmin=", rightmin
    scale = gPad.GetUymax()/rightmax


    for i in range(0, grM.GetN()):
        grM.GetY()[i]*=scale


    grM.Draw("same")



    axis = TGaxis(gPad.GetUxmax(),gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(), 0, rightmax, 510,"+L")
    axis.SetLineColor(kGreen+2)
    axis.SetTextColor(kGreen+2)
    axis.SetLabelColor(kGreen+2)
    axis.SetTitleOffset(1.5)

    axis.Draw()
    axis.SetTitle("Mobility [cm^{2}/Vs]")
    gPad.RedrawAxis()
    gPad.Update()
    gPad.Modified()

    
    bla=raw_input()

