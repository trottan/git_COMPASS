import ROOT
import numpy as np
import math


fname = "merged*.root"
rdf = ROOT.RDataFrame("tree",fname)



rdf = rdf.Filter("Q2 > 1 && Q2 < 10")
rdf = rdf.Filter("W > 5")
rdf = rdf.Filter("y > 0.1 && y < 0.9")
rdf = rdf.Filter("pt2 < 0.5 && pt2 > 0.01")
rdf = rdf.Filter("M_rho0 > 0.5 && M_rho0 < 1.1")
rdf = rdf.Filter("mom_mv > 15.0")



rdf = rdf.Define("vals","""
        auto Mp = 0.938272;
        auto Mpi = 0.13957;
        auto Mmu = 0.000105;
        auto Pmu_in = TMath::Sqrt(mu_in[0]*mu_in[0]+mu_in[1]*mu_in[1]+mu_in[2]*mu_in[2]);
        auto Mmu_in = TMath::Sqrt(mu_in[3]*mu_in[3]-Pmu_in*Pmu_in);

        auto Pmu_out = TMath::Sqrt(mu_out[0]*mu_out[0]+mu_out[1]*mu_out[1]+mu_out[2]*mu_out[2]);
        auto Mmu_out = TMath::Sqrt(mu_out[3]*mu_out[3]-Pmu_out*Pmu_out);


        TLorentzVector targ(0,0,0,Mp); //target proton

        TLorentzVector pip,pim;
        pip.SetPxPyPzE(pi_1[0],pi_1[1],pi_1[2],pi_1[3]);
        pim.SetPxPyPzE(pi_2[0],pi_2[1],pi_2[2],pi_2[3]);
        TLorentzVector beam, ele;
        beam.SetPxPyPzE(mu_in[0],mu_in[1],mu_in[2],mu_in[3]);
        ele.SetPxPyPzE(mu_out[0],mu_out[1],mu_out[2],mu_out[3]);

        auto q = beam - ele;
        auto rho = pip +pim;
        auto pro = beam+targ-pip-pim-ele;
auto Rot_Matrix = [&](TLorentzVector vector, int Lab2CM_or_CM2Lab, double Theta_Rot, double Phi_Rot){ 

        double Rot_X1 = vector.X();
        double Rot_Y1 = vector.Y();
        double Rot_Z1 = vector.Z();

        double Rot_X = Rot_X1;
        double Rot_Y = Rot_Y1;
        double Rot_Z = Rot_Z1;


        // Lab2CM_or_CM2Lab is a parameter which determines if you rotating from the lab frame to the CM frame, or if you are rotating back in the opposite direction
        // Lab2CM_or_CM2Lab = -1 gives a rotation to the CM frame (from the lab frame)
        // Lab2CM_or_CM2Lab = +1 gives a rotation to the lab frame (from the CM frame)


        Theta_Rot = -1*Theta_Rot;   // Always give the angle of rotation Theta as the value given by .Theta()
                                    // This subroutine will handle the fact that the matrix rotation wants the negative of the angle of rotation


        // Rotation to Lab Frame
        if(Lab2CM_or_CM2Lab == -1){
            Rot_X = Rot_X1*TMath::Cos(Theta_Rot)*TMath::Cos(Phi_Rot) - Rot_Z1*TMath::Sin(Theta_Rot) + Rot_Y1*TMath::Cos(Theta_Rot)*TMath::Sin(Phi_Rot);
            Rot_Y = Rot_Y1*TMath::Cos(Phi_Rot) - Rot_X1*TMath::Sin(Phi_Rot);
            Rot_Z = Rot_Z1*TMath::Cos(Theta_Rot) + Rot_X1*TMath::Cos(Phi_Rot)*TMath::Sin(Theta_Rot) + Rot_Y1*TMath::Sin(Theta_Rot)*TMath::Sin(Phi_Rot);
        }


        // Rotation to CM Frame
        if(Lab2CM_or_CM2Lab == 1){
            Rot_X = Rot_X1*TMath::Cos(Theta_Rot)*TMath::Cos(Phi_Rot) + Rot_Z1*TMath::Cos(Phi_Rot)*TMath::Sin(Theta_Rot) - Rot_Y1*TMath::Sin(Phi_Rot);
            Rot_Y = Rot_Y1*TMath::Cos(Phi_Rot) + Rot_X1*TMath::Sin(Phi_Rot)*TMath::Cos(Theta_Rot) + Rot_Z1*TMath::Sin(Theta_Rot)*TMath::Sin(Phi_Rot);
            Rot_Z = Rot_Z1*TMath::Cos(Theta_Rot) - Rot_X1*TMath::Sin(Theta_Rot);
        }



        TLorentzVector vector_Rotated(Rot_X, Rot_Y, Rot_Z, vector.E());

        return vector_Rotated;


    };

    double Theta_q = q.Theta();
    double Phi_el = ele.Phi();

    // CM rotation needed

    auto beam_Clone = Rot_Matrix(beam, -1, Theta_q, Phi_el);
    auto ele_Clone  = Rot_Matrix(ele,  -1, Theta_q, Phi_el);
    auto pip_Clone = Rot_Matrix(pip,  -1, Theta_q, Phi_el);
    auto q_Clone = Rot_Matrix(q,  -1, Theta_q, Phi_el);
    auto rho_Clone = Rot_Matrix(rho,  -1, Theta_q, Phi_el);
    auto pro_Clone = Rot_Matrix(pro,  -1, Theta_q, Phi_el);
    auto targ_Clone = Rot_Matrix(targ, -1, Theta_q, Phi_el);
    
    auto fCM = q_Clone + targ_Clone;
    auto boost = -(fCM.BoostVector());
    //TVector3 boost(-fCM.Px()/fCM.E(),-fCM.Py()/fCM.E(),-fCM.Pz()/fCM.E());


    auto qlvBoost(q_Clone);
    auto eleBoost(ele_Clone);
    auto beamBoost(beam_Clone);
    auto targBoost(targ_Clone);
    auto pipBoost(pip_Clone);
    auto rhoBoost(rho_Clone);
    auto proBoost(pro_Clone);

    qlvBoost.Boost(boost);
    eleBoost.Boost(boost);
    beamBoost.Boost(boost);
    targBoost.Boost(boost);
    pipBoost.Boost(boost);
    rhoBoost.Boost(boost);
    proBoost.Boost(boost);
    
    //Get Big Phi -Angle between lepton scattering plane and rho prodution plane in CM frame
    auto A = (qlvBoost.Vect().Cross(rhoBoost.Vect())).Dot((beamBoost.Vect().Cross(eleBoost.Vect())));
    auto B = (qlvBoost.Vect().Cross(rhoBoost.Vect())).Mag()*(beamBoost.Vect().Cross(eleBoost.Vect())).Mag();
    auto c0 = (qlvBoost.Vect().Cross(eleBoost.Vect())).Dot(rhoBoost.Vect());
    auto Phi =(c0/TMath::Abs(c0)) * TMath::ACos(A/B);//*TMath::RadToDeg();
    

    //Rho rest frame rotation Need

    auto Phi_rho =  rhoBoost.Phi();
    auto Theta_rho = rhoBoost.Theta();

    auto pip_Clone_rf = Rot_Matrix(pipBoost,  -1, Theta_rho, Phi_rho);
    auto rho_Clone_rf = Rot_Matrix(rhoBoost,  -1, Theta_rho, Phi_rho);
    auto q_Clone_rf = Rot_Matrix(qlvBoost,  -1, Theta_rho, Phi_rho);
    auto pro_Clone_rf = Rot_Matrix(proBoost,  -1, Theta_rho, Phi_rho);


    auto fRF =  -rho_Clone_rf;
    auto boost2 = -(fRF.BoostVector());

    auto pipBoost2(pip_Clone_rf);
    auto rhoBoost2(rho_Clone_rf);
    auto qBoost2(q_Clone_rf);
    auto proBoost2(pro_Clone_rf);




    pipBoost2.Boost(boost2);
    rhoBoost2.Boost(boost2);
    qBoost2.Boost(boost2);
    proBoost2.Boost(boost2);


    //Get Small phi -Azimuthal angle between rho production plane and decay plan in vector-meson rest frame
    A =  (qlvBoost.Vect().Cross(rhoBoost.Vect())).Dot((rhoBoost.Vect().Cross(pipBoost.Vect())));
    B = (qlvBoost.Vect().Cross(rhoBoost.Vect())).Mag() * (rhoBoost.Vect().Cross(pipBoost.Vect())).Mag();
    auto a2 = (qlvBoost.Vect().Cross(rhoBoost.Vect())).Dot(pipBoost.Vect());
    auto phi = (c0/TMath::Abs(c0)) *TMath::ACos(A/B);//*TMath::RadToDeg();
    auto Bphi = B;





    //Get Cos(theta) - Decay of Pi+ in the vector meson rest Frame
    A = proBoost2.Vect().Dot(pipBoost2.Vect());
    B = proBoost2.Vect().Mag()*pipBoost2.Vect().Mag();
    auto cosTheta = -A/B;

    double gamma = 2.0*Mp*Xbj/TMath::Sqrt(Q2);
    double epsA = 1.0-y-1.0/4.0*gamma*gamma*y*y;
    double epsB = 1.0-y+1.0/2.0*y*y+1.0/4.0*gamma*gamma*y*y;
    double eps = epsA/epsB;

    return vector<double> {Phi,phi,cosTheta,eps,Bphi,a2};




""")

rdf = rdf.Define("LzVec_Print","""

auto TempPart = LzVecOut;
return vector<double>  {TempPart.Px(),TempPart.Py(),TempPart.Pz(),TempPart.E(),TempPart.X(),TempPart.Y(),TempPart.Z(),TempPart.M()}
""")

values = ["vals"]
#values = ["mu_in"]
#output = rdf.Display(values,10).GetValue().AsString()
#output = rdf.Display(values,  rdf.Count().GetValue()).GetValue().AsString()

#print(output)
rdf = rdf.Define("Phi","vals[0]")
rdf = rdf.Define("phi","vals[1]")
rdf = rdf.Define("cosTheta","vals[2]")
rdf = rdf.Define("eps","vals[3]")

def makeRoot(rdf):
    output_file = "Kin_eppippim_MC.root"
    output = ROOT.TFile(output_file, "RECREATE")
    columns_to_keep = ["Phi","phi","cosTheta","eps"]  # Only keep these
    rdf.Snapshot("tree", output_file, columns_to_keep)
    # Close the output file
    output.Close()

#makeRoot(rdf)






def makePlots(rdf):
    hPhi = rdf.Histo1D(("hPhi", "#Phi;#Phi",200,-4,4), "Phi")
    hphi = rdf.Histo1D(("hphi", "#phi;#phi",200,-4,4), "phi")
    hcosTheta = rdf.Histo1D(("hcosTheta", "cos(#theta);cos(#theta)",200,-1.5,1.5), "cosTheta")



    c1 = ROOT.TCanvas("c1","c1", 1000,700)
    #c1.Divide(2,2,.01,0.01)
    c1.Draw()
    c1.Print("Angles.pdf[")
    #c1.cd(1)
    hPhi.Draw()
    c1.Print("Angles.pdf")
    #c1.cd(2)
    hphi.Draw()
    c1.Print("Angles.pdf")
    #c1.cd(3)
    hcosTheta.Draw()
    c1.Print("Angles.pdf")
    c1.SaveAs("Angles.pdf]")
makePlots(rdf)
