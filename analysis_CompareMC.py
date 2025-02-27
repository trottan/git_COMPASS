#!/usr/bin/env python3
import ROOT
fname_MC = "rho_2016_MC/merged*.root"
rdf_MC = ROOT.RDataFrame("tree",fname_MC)

Q2Cut = "Q2 > 1 && Q2 < 10"
WCut = "W > 5"
yCut = "y > 0.1 && y < 0.9"
PtCut = "pt2 < 0.5 && pt2 > 0.01"
mrhoCut = "M_rho0 > 0.5 && M_rho0 < 1.1"
prhoCut = "mom_mv > 15.0"

cutW_MC = rdf_MC.Filter(Q2Cut).Filter(yCut).Filter(PtCut).Filter(mrhoCut).Filter(prhoCut)
cutM_MC = rdf_MC.Filter(Q2Cut).Filter(yCut).Filter(PtCut).Filter(WCut).Filter(prhoCut)
cutPt_MC = rdf_MC.Filter(Q2Cut).Filter(yCut).Filter(mrhoCut).Filter(WCut).Filter(prhoCut)
cutQ2_MC = rdf_MC.Filter(PtCut).Filter(yCut).Filter(PtCut).Filter(WCut).Filter(prhoCut)
cutAll_MC = rdf_MC.Filter(PtCut).Filter(yCut).Filter(PtCut).Filter(WCut).Filter(prhoCut).Filter(Q2Cut)

fname = "rho_2016/merged*.root"
rdf = ROOT.RDataFrame("tree",fname)


cutW = rdf.Filter(Q2Cut).Filter(yCut).Filter(PtCut).Filter(mrhoCut).Filter(prhoCut)
cutM = rdf.Filter(Q2Cut).Filter(yCut).Filter(PtCut).Filter(WCut).Filter(prhoCut)
cutPt = rdf.Filter(Q2Cut).Filter(yCut).Filter(mrhoCut).Filter(WCut).Filter(prhoCut)
cutQ2 = rdf.Filter(PtCut).Filter(yCut).Filter(PtCut).Filter(WCut).Filter(prhoCut)
cutAll  = rdf.Filter(PtCut).Filter(yCut).Filter(PtCut).Filter(WCut).Filter(prhoCut).Filter(Q2Cut)








def makePlots(rdf,rdf_MC,outputName):
    #ans = input("Would you like plots?(Q2, MM and MGamma)")
    #enable = True if (ans == "yes" or ans == "y") else False
    if(True):
        hists = []
        hists_MC = []
        hq2 = rdf.Histo1D(("hq2", "Q^{2};Q^{2} [GeV^{2}]",400,0,10), "Q2")
        hMRho = rdf.Histo1D(("hMrho", "M(#pi+#pi-);Invariant Mass #pi+#pi- [GeV]",400,0,1.5), "M_rho0")
        hq2xb = rdf.Histo2D(("hq2xb", "Q^{2} vs x_{B};x_{B};Q^{2}", 100,0,0.1,100,0,11), "Xbj", "Q2")
        hww = rdf.Histo1D(("hW", "W;W [GeV^{2}]",400,0,20), "W")
        hpt = rdf.Histo1D(("hpt", "p_{t}^2;p_{t}2 [GeV]",400,0,2), "pt2")
        hEmiss = rdf.Histo1D(("hEmiss", "EMiss;EMiss [GeV]",400,-10,40), "Emiss")





        hq2_MC = rdf_MC.Histo1D(("hq2", "Q^{2};Q^{2} [GeV^{2}]",400,0,10), "Q2","weight_all")
        hMRho_MC = rdf_MC.Histo1D(("hMrho", "M(#pi+#pi-);Invariant Mass #pi+#pi- [GeV]",400,0,1.5), "M_rho0","weight_all")
        hq2xb_MC = rdf_MC.Histo2D(("hq2xb", "Q^{2} vs x_{B};x_{B};Q^{2}", 100,0,0.1,100,0,11), "Xbj", "Q2","weight_all")
        hww_MC = rdf_MC.Histo1D(("hW", "W;W [GeV^{2}]",400,0,20), "W","weight_all")
        hpt_MC = rdf_MC.Histo1D(("hpt", "p_{t}^2;p_{t}2 [GeV]",400,0,2), "pt2","weight_all")
        hEmiss_MC = rdf_MC.Histo1D(("hEmiss", "EMiss;EMiss [GeV]",400,-10,40), "Emiss","weight_all")
        #hMM = rdf.Histo1D(("hMM","Missing Mass^2",100, -20,20),"MM")
        #hMGammaGamma = rdf.Histo1D(("hMGammaGamma","\gamma\gamma Mass",100, -0.6,0.6),"mGammaGamma")
        #fit1d(hMGammaGamma)
        hists.append(hq2)
        hists.append(hMRho)
        hists.append(hww)
        hists.append(hpt)
        hists.append(hEmiss)
        #hists.append(hq2xb)

        hists_MC.append(hq2_MC)
        hists_MC.append(hMRho_MC)
        hists_MC.append(hww_MC)
        hists_MC.append(hpt_MC)
        hists_MC.append(hEmiss_MC)
        #hists_MC.append(hq2xb_MC)

        c1 = ROOT.TCanvas("c1","c1", 1200,900)
        #c1.Divide(3,2,0.0001,0.0001)
        c1.SetLogx()
        c1.Draw()
        c1.Print(outputName + '[')
        for i,j in zip(hists,hists_MC):
            #if isinstance(i, ROOT.TH1) and not isinstance(i, ROOT.TH2):
            #    print("This is a TH1 (1D histogram).")
            #    i.Draw()
            #elif isinstance(i, ROOT.TH2):
            #    print("This is a TH2 (2D histogram).")
            #    i.Draw("colz")
                
            maxNum = max(rdf.Count().GetValue(),rdf_MC.Count().GetValue())
            i.GetYaxis().SetRangeUser(0,500)
            i.DrawNormalized()
            j.SetLineColor(2)
            j.DrawNormalized("same")
            c1.Print(outputName)
        c1.Print(outputName + ']')


#makePlots(cutAll,cutAll_MC,"compareMC.pdf")
#makePlots(cutW,cutW_MC,"compareMC_W.pdf")
#makePlots(cutM,cutM_MC,"compareMC_M.pdf")
#makePlots(cutPt,cutPt_MC,"compareMC_Pt.pdf")
makePlots(cutQ2,cutQ2_MC,"compareMC_Q2.pdf")
