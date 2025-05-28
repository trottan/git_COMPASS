#include "iostream"
#include "cmath"
#include "sstream"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "Phast.h"
#include "PaEvent.h"
#include "PaAlgo.h"
#include "PaVertex.h"
#include "G3part.h"
#include "PaParticle.h"
#include "PaTrack.h"
#include "PaTPar.h"
#include "PaSetup.h"
#include "TROOT.h"
#include "PaMChit.h"
#include <boost/foreach.hpp>
#include "string.h"
#include "TargetCell.h"
#include "Toololols_KinFitterInterface.hh"
#include "Tools_Camera.hh"
#include "TiS_range.cc"
#include "UConn_Tools.h"


bool verbose_mode = true; // Create an instance of verbose_mode
static PaCamera* cam_inst;
PaHodoHelper* HodoHelper; 
EventFlags eventFlags;

void UserEvent610(PaEvent& e)
{
    //  NLUDATA ld;                   // create NLUDATA structure
    //  vector<LUJET> lujets;         // create LUJET structure
    //  int nparticles;
    //  e.MCgen(ld);
    //  e.MCgen(nparticles,lujets);

  Bool_t intarget;
  static PaHodoHelper* HodoHelper;
  //static TiSRange*     tis_range;
  const double M_mu = G3partMass[5];  // muon   mass
  const double M_pi = G3partMass[8];  // pion   mass
  const double M_p = G3partMass[14];  // proton mass

  static TTree* tree(NULL); // tree for booking

  long long     UniqueEvNum;                          // number of event
  static int    Run, Spill, Year, EvNum;              // event properties
  static int    LastRun = -1;
  static int    trig;                                 // trig mask
  static double vx, vy, vz;                           // X coordinate of primary vertex
  static double Xprim, Yprim, Zprim;                  // X coordinate of primary vertex
  static int    Nout, NVertex, NGoodCandidates;       // Number of tracks in primary vertex, number of vertices per event, number of good candiates in CAMERA
  static int    trig_mask;
  static double Chi2Vtxt, Chi2Vtxt_ndf;               // Chi2 of primary vertex
  static double y, W, Q2, Xbj, pt2;                   // kinematic variables
  static float  tbeam;                                // time of the incoming beam
  static int    qbeam;                                // charge of the incoming beam
  static float  qq_pions;                             // (h1.Q() + h2.Q())
  static float  mom_mv;                               // meson vector momentum
  static float  M2pi, epsm, polbm;
  double        nu, E, E_scat;                        // nu, incoming muon energy, scattered muon energy
  double        mom, mom_err;                         // incoming muon momentum
  int           charge_in, charge_out;                // incoming muon charge, outgoing muon charge
  int           cl1_EC0, cl1_EC1, cl1_EC2;            // first cluster in ECALs
  int           cl2_EC0, cl2_EC1, cl2_EC2;            // first cluster in ECALs
  double        M_rho0, E_rho0, cl1_Energy, cl2_Energy;
  double        cl1_sigmaT, cl2_sigmaT, cl1_time, cl2_time, cl1_x, cl2_x, cl1_y, cl2_y, cl1_z, cl2_z;
  double        Emiss, Mx2;                           // kinematic variables
  static float  OpAngle;                              // opening angle
//  double        weight_all; 
   static double time_in_spill;
   static double MM;
  Double_t mu_in[4];
  Double_t mu_out[4];
  Double_t pi_1[4];
  Double_t pi_2[4];
  Double_t pro[4];

  TLorentzVector LzVecMu0;   // beam mu Lorentz vector
  TLorentzVector LzVecMu1;   // scattered mu Lorentz vector
  TLorentzVector LzVecOut1;  // pi+
  TLorentzVector LzVecOut2;  // pi-
  TLorentzVector LzVecOut;   // pi+ + pi-
//  TLorentzVector LzVecPi0;   // pi0
  TLorentzVector LzVecPspec;
  TLorentzVector LzVec_p_rec[200];
  TLorentzVector LzVecPro;
  


static double phi_MM;
static double phi_cam;
static double pt_MM;
static double pt_cam;
static double delta_Z;
static double dPhi;
static double dpt;
static double vertex_x;
static double vertex_y;
static double vertex_z;
TVector3 R_vertex;
TVector3 pos_in_ring_b[200];
TVector3 pos_in_ring_a[200];
 // int Nout_gen;i
   static TargetCell& target = TargetCell::Instance();
   bool TiS_flag  = false;
   bool trig_flag = false;
   bool flux_flag  = false;
   bool outMu_flag = false;  
   bool mu_flag = false;  
   bool targ_flag = false;  
   eventFlags.createFlag("Only 3 particles", "Tracks only have three outgoing particles");
   eventFlags.createFlag("Pen_Len_had", "Penetration length of Hadron");
   eventFlags.createFlag("Qu_Fit_Had", "Hadrons have good Quality of fit");
   eventFlags.createFlag("HAD_SM1", "The track of the hadron is before the first magnet");
   eventFlags.createFlag("Opp_charge", "Hadrons have opposite charge");
   eventFlags.createFlag("Emiss_flag", "Missing energy cut");
   eventFlags.createFlag("EE_flag", "Scattered Muon Energy less than beam");
   eventFlags.createFlag("W_flag", "W cut");
   eventFlags.createFlag("Q2_flag", "Q2 cut");
   eventFlags.createFlag("pt_flag", "pt2 cut");
   eventFlags.createFlag("nu_flag", "nu cut");
   eventFlags.createFlag("rhoP_flag", "rho momentum cut");
   eventFlags.createFlag("y_flag", "y cut");
   eventFlags.createFlag("Minv_flag", "Invarnt Mass cut");
   eventFlags.resetFlags(); // Reset all event statistic counter flags to false
////////////////////////////////////HISTOGRAM BOOKING/////////////////////////////////////////////////////

  static bool first(true);
  if(first){ // histograms booking block

	 HodoHelper=&PaHodoHelper::Init("/afs/cern.ch/user/n/ntrotta/workDir/work_phast_DVCS2016/phast/trigger_config/2016"); 
    tree = new TTree("tree","Variables");
    tree->Branch("UniqueEvNum",         &UniqueEvNum,                   "UniqueEvNum/L");
    tree->Branch("Year",                  &Year,                        "Year/I");
    tree->Branch("Run",                           &Run,                           "Run/I");
    tree->Branch("Spill",                       &Spill,                         "Spill/I");
    tree->Branch("EvNum",                 &EvNum,                       "Evnum/I");
    tree->Branch("trig",                          &trig,                          "trig/I");
    tree->Branch("qbeam",               &qbeam,               "qbeam/I");
    tree->Branch("intarget",                    &intarget,                      "intarget/O");
    tree->Branch("vx",                          &vx,                          "vx/D");
    tree->Branch("vy",                          &vy,                                "vy/D");
    tree->Branch("vz",                    &vz,                        "vz/D");
    tree->Branch("Xprim",                         &Xprim,                           "Xprim/D");
    tree->Branch("Yprim",                         &Yprim,                                 "Yprim/D");
    tree->Branch("Zprim",                         &Zprim,                                 "Zprim/D");
    tree->Branch("Chi2Vtxt",                    &Chi2Vtxt,                      "Chi2Vtxt/D");
    tree->Branch("Chi2Vtxt_ndf",              &Chi2Vtxt_ndf,                "Chi2Vtxt_ndf/D");
    tree->Branch("NVertex",                 &NVertex,                         "NVertex/I");
    tree->Branch("Nout",                        &Nout,                            "Nout/I");
 //   tree->Branch("Nout_gen",                  &Nout_gen,                "Nout_gen/I");
    tree->Branch("y",                           &y,                               "y/D");
    tree->Branch("W",                           &W,                             "W/D");
    tree->Branch("Q2",                          &Q2,                              "Q2/D");
    tree->Branch("Xbj",                         &Xbj,                             "Xbj/D");
    tree->Branch("pt2",                         &pt2,                             "pt2/D");
    tree->Branch("tbeam",                 &tbeam,                         "tbeam/F");
    tree->Branch("qq_pions",              &qq_pions,                    "qq_pions/F");
    tree->Branch("nu",                    &nu,                            "nu/D");
    tree->Branch("E",                           &E,                               "E/D");
    tree->Branch("E_scat",              &E_scat,                      "E_scat/D");
    tree->Branch("mom",                   &mom,                         "mom/D");
    tree->Branch("mom_err",             &mom_err,                       "mom_err/D");
    tree->Branch("epsm"    ,              &epsm,                        "epsm/F");
    tree->Branch("polbm",                       &polbm,                           "polbm/F");
    tree->Branch("M2pi",                        &M2pi,                            "M2pi/F");
    tree->Branch("mom_mv",                      &mom_mv,                        "mom_mv/F");
    tree->Branch("charge_in",               &charge_in,             "charge_in/I");
    tree->Branch("charge_out",              &charge_out,          "charge_out/I");
    tree->Branch("Mx2",                           &Mx2,                 "Mx2/D");
//    tree->Branch("M_pi0",                             &M_pi0,                     "M_pi0/D");
    tree->Branch("M_rho0",                            &M_rho0,                  "M_rho0/D");
    tree->Branch("E_rho0",                      &E_rho0,                        "E_rho0/D");
    tree->Branch("NGoodCandidates",     &NGoodCandidates,       "NGoodCandidates/I");
    tree->Branch("Emiss",                         &Emiss,                           "Emiss/D");
    tree->Branch("OpAngle",                     &OpAngle,                       "OpAngle/F");
    tree->Branch("mu_in",               mu_in,                "mu_in[4]/D");
    tree->Branch("mu_out",              mu_out,               "mu_out[4]/D");
    tree->Branch("pi_1",                pi_1,                 "pi_1[4]/D");
    tree->Branch("pi_2",                pi_2,                 "pi_2[4]/D");
   // tree->Branch("pro",                pro,                 "pro[4]/D");
    tree->Branch("LzVecMu0",                  &LzVecMu0);
   // tree->Branch("LzVecPi0",                &LzVecPi0);
    tree->Branch("LzVecMu1",              &LzVecMu1);
    tree->Branch("LzVecOut",                  &LzVecOut);
    tree->Branch("LzVecOut1",                 &LzVecOut1);
    tree->Branch("LzVecOut2",                 &LzVecOut2);
    tree->Branch("LzVecPspec",                &LzVecPspec);
  //  tree->Branch("weight_all",                      &weight_all,              "weight_all/D");
    first=false;

   target.Init(e);
  } // end of histogram booking



/////////////////////////////////////END OF HISTOGRAM BOOKING///////////////////////////////////////////////////////

/*enum trigger {   
                Tiger = 1<<0,
                MT = 1<<1,
                LT = 1<<2,
                OT = 1<<3,
                CT = 1<<4,
                IV = 1<<5,
                HaloT = 1<<6,
                BT = 1<<7,
                Tiger_only = 1<<8,
                LAST = 1<<9,
                TRand = 1<<10,
                NRand = 1<<11
    };

    trig_mask = e.TrigMask();
    trig_mask = trig_mask&2047;

    std::string trigCheck = ""; // empty string to store trigger information for debugging 
    if (trig_mask & MT) {
        trigCheck += "MT ";
        trig_flag = true;
    }
    if (trig_mask & LT) {
        trigCheck += "LT ";
        trig_flag = true;
    }
    if (trig_mask & OT) {
        trigCheck += "OT ";
        trig_flag = true;
    }
    if (trig_mask & LAST) {
        trigCheck += "LAST ";
        trig_flag = true;
    }*/
//****************************************************************************************************************************
//***************************************************** BEGGINING OF CUTS *****************************************************
//****************************************************************************************************************************

  PaVertex tempVertex; // vertex radial position
  PaTrack tempTrack;   // vertex longitudinal position

  Year = e.Year(); // event properties
  Run = e.RunNum();
  Spill = e.SpillNum();
  EvNum = e.EvInSpill();
  UniqueEvNum = e.UniqueEvNum();

  cam_inst=&PaCamera::GetInstance();
  cam_inst->NewEvent(e);//,e.IsMC());
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
time_in_spill= e.TimeInSpill();
if (Run != LastRun) { // Reinitialize HodoHelper and tis_range only if the run number changes 
       // HodoHelper = & PaHodoHelper::Init("", true);  
        //tis_range  = new TiSRange("/afs/cern.ch/user/n/ntrotta/workDir/work_phast_DVCS2016/phast/user/flux_files/flux_Johannes/2016/flux_files");
	Set_TiSrange("/afs/cern.ch/user/n/ntrotta/workDir/work_phast_DVCS2016/phast/user/flux_files/flux_Johannes/2016/flux_files",Run,Run);
	LastRun = Run;  // Update LastRun to the current run number	
}
eventFlags.setFlagByName("allEvts_flag", true);
TiS_flag = Check_TiS_window(Run,Spill,time_in_spill);



for (int iv = 0; iv < e.NVertex(); iv++) { // begin loop over vertices
		
	const PaVertex & v = e.vVertex(iv);
	Nout = v.NOutParticles(); // number of tracks in vertex	
	
	if (!v.IsPrimary()) continue;
	eventFlags.setFlagByName("pVtx_flag", true);
	if (!v.IsBestPrimary()) continue;	
	if(Nout !=3) continue;
	eventFlags.setFlagByName("Only 3 particles", true);

	Zprim = v.Pos(2);
	Yprim = v.Pos(1);
	Xprim = v.Pos(0);
	
	static PaParticle beam; 
      	static PaTrack beam_track; 
      	static PaTPar Par_beam;
	static BeamFluxParams beamParams; // Create an instance of BeamFluxParams
	//flux_flag = beamFluxCheck(e, v, iv, Run, true, beamParams, beam, beam_track, Par_beam, eventFlags);
	//if (!flux_flag) continue;
	static PaParticle outMu;
        static PaTrack outMu_track;
        static PaTPar Par_outMu;
        
	static OutMuParams outMuParams; // Create an instance of OutMuParams

     mu_flag = crossCheck2(e,v,iv,Run,beamParams, beam, beam_track, Par_beam, HodoHelper,outMuParams, outMu, outMu_track, Par_outMu,eventFlags);
       if(!mu_flag) {continue;}

        PaVertex tempVertex = v;
        PaTrack tempTrack = beam_track;
	targ_flag = true;
        /*if(!(target.InTarget(&tempVertex, 1.9))) {
                targ_flag = false;
        }
        if(!target.InTarget(&tempTrack,-318.5,-78.5,1.9)){
                targ_flag = false;
        }*/
       	/*if(!(targ_flag)){continue;}
	eventFlags.setFlagByName("vtxInTarget_flag", true);
	*/
	
	qbeam= beam.Q();
	mom = beam_track.vTPar(0).Mom();	
	
	/*if (!( ((e.TrigMask() & 0x2) != 0) || ((e.TrigMask() & 0x4) != 0) || ((e.TrigMask() & 0x8) != 0) || ((e.TrigMask() & 0x200) != 0) )){
	       //trig_flag = true;
	       continue;
	}*/
	
	
	
	/*static PaParticle outMu; 
      	static PaTrack outMu_track; 
      	static PaTPar Par_outMu;
      	static OutMuParams outMuParams; // Create an instance of OutMuParams*/
      	//outMu_flag = outMuCheck(e, v, iv, Run, beam, HodoHelper, true, outMuParams, outMu, outMu_track, Par_outMu, eventFlags);
	//if(!outMu_flag){continue;}
	
	int imu1 = HodoHelper->iMuPrim(v,false,false,true,true,15.0,true,true); //v.iMuPrim(false,true,true,false,15);
	
		
	LzVecMu0 = Par_beam.LzVec(M_mu); // calculate beam mu Lorentz vector
     	LzVecMu1 = Par_outMu.LzVec(M_mu); // calculate scattered mu Lorentz vector  
        
     	TVector3 VecMu0 = LzVecMu0.Vect();
      	TVector3 VecMu1 = LzVecMu1.Vect();
      	TVector3 scattering_plane = VecMu0.Cross(VecMu1);
      
      	y   = (LzVecMu0.E()-LzVecMu1.E())/LzVecMu0.E(); // fraction of lepton energy lost in the laboratory system  
      	Q2  = PaAlgo::Q2 (LzVecMu0, LzVecMu1); // virtuality of the photon
      	Xbj = PaAlgo::xbj(LzVecMu0, LzVecMu1); // Bjorken scaling variable
      	W   = sqrt(PaAlgo::W2 (LzVecMu0, LzVecMu1)); // invariant mass of the photon-nucleon system
      
      	tbeam = beam_track.MeanTime(); // time of the incoming beam
     	/*if (Q2 < 0.5) continue;  
	eventFlags.setFlagByName("Q2_DIS_flag", true);*/
	//if(Nout !=3) continue;
	//eventFlags.setFlagByName("Only 3 particles", true);

//****************************** Outging particles ******************************************************************************** 

      TLorentzVector LzVecOut1(0,0,0,0); // for  4-vector  of hadron 1   
      TLorentzVector LzVecOut2(0,0,0,0); // for  4-vector  of hadron 2

      int kmu=-1;
      int hd1=-1;
      int hd2=-1;
      for(int j = 0; j < v.NOutParticles(); j++) { //loop over outgoing paticles of the vertex 
      int ip = v.iOutParticle(j);  // index of outgoing particle

      if(ip==imu1) {kmu=j;
        continue;}      // skip mu'             
      if(hd1 == -1){ hd1=j;
        continue;}               
      if(hd2 == -1){ hd2=j;                
        continue;}           
      }           

       int ip1 = v.iOutParticle(hd1);
       int ip2 = v.iOutParticle(hd2);
       PaParticle h1 = e.vParticle(ip1);
       PaParticle h2 = e.vParticle(ip2);

       const PaTrack & Trackp1 = e.vTrack(h1.iTrack());
       const PaTrack & Trackp2 = e.vTrack(h2.iTrack());
    if(Trackp1.ZLast()<350 || Trackp2.ZLast()<350) continue; 
    if(Trackp1.XX0() > 10) continue; // penetration length of hadron track shorter than 10 radiation lengths
    if(Trackp2.XX0() > 10) continue; // penetration length of hadron track shorter than 10 radiation lengths
	eventFlags.setFlagByName("Pen_Len_had", true);
    if(Trackp1.Chi2tot()/Trackp1.Ndf() > 10) continue; // good fit quality of outgoing hadron ch_red^2 < 10
    if(Trackp2.Chi2tot()/Trackp2.Ndf() > 10) continue; // good fit quality of outgoing hadron ch_red^2 < 10
       eventFlags.setFlagByName("Qu_Fit_Had", true);
    if(Trackp1.ZFirst()>350 || Trackp2.ZFirst()>350) continue;// tracks of hadrons start upstream from SM1
     eventFlags.setFlagByName("HAD_SM1", true);

    //For rho hadrons of oppsite charge
    /*if(h1.Q() == h2.Q()) continue; //skip equaly charged particles 
	eventFlags.setFlagByName("Opp_charge", true);
    if(h1.Q() < h2.Q()) swap(h1,h2); //now first particle is positive            
    */   
     
    // For lepto Comparison

     if(h1.Q() != h2.Q()) continue;
       
       qq_pions= h1.Q()+h2.Q() ;      
       LzVecOut1 =   h1.ParInVtx(iv).LzVec(M_pi);  //  4-vector  of hadron 1 positive
       LzVecOut2 =   h2.ParInVtx(iv).LzVec(M_pi);  //  4-vector  of hadron 2 positive             

       polbm = PaAlgo::GetBeamPol(mom, 2011);
       
       E= LzVecMu0.E(); //beam muon energy
       E_scat = LzVecMu1.E(); //scattered muon energy
       nu = E - E_scat;  
       epsm =(1-y-y*y*Q2/(4.0*nu*nu))/(1-y+0.25*y*y*(Q2/(nu*nu) +2.0));

//****************************** Proton **************************************************************************************************
	  M2pi = (LzVecOut1+LzVecOut2).Mag();
	  TLorentzVector q =  LzVecMu0 - LzVecMu1; // virtual photon's 4-vector
	  TLorentzVector LzVecP0(0,0,0,M_p); //target proton
	  TLorentzVector LzVecOut = LzVecOut1 + LzVecOut2;
	  TLorentzVector LzVecPspec= LzVecMu0 + LzVecP0 - LzVecMu1 - LzVecOut; //predicted recoil proton from spectrometer
	  mom_mv= (LzVecOut.Vect()).Mag();
	  double Mmiss2 = LzVecPspec.Mag2();
	  Emiss = (Mmiss2-(M_p*M_p))/(2*M_p);
	  M_rho0 = LzVecOut.M();
	  E_rho0 = LzVecOut.E();
	//  M_pi0 = LzVecPi0.M();
	  TVector3 qq3=q.Vect();
    	  pt2=LzVecOut.Perp2(qq3);

	  
//Proton using Camera (Code is still in progress) 

  // Proton from Camera	
 	/*	vector <CameraProton> protons=cam_inst->GetGoodCandidates(v);
		R_vertex.SetXYZ(v.Pos(0),v.Pos(1),v.Pos(2));
		bool tot_flag = false;

		vertex_x= R_vertex.X();
		vertex_y= R_vertex.Y();
		vertex_z= R_vertex.Z(); 
		int count_pro =0;
		int num_comb = 0;
		MM = -9999;
		delta_Z = -9999;
		dPhi = -9999;
		dpt = -9999;
		
		for (unsigned int i=0;i<protons.size();i++){
			double beta;
			beta=protons[i].beta;
			if(beta < 0.1 || beta > 1) continue;
			LzVec_p_rec[num_comb]=protons[i].p4;
			if(LzVec_p_rec[num_comb].Mag() ==0) continue;
			
			pos_in_ring_a[num_comb]=protons[i].Ahit.vec;
			pos_in_ring_b[num_comb]=protons[i].Bhit.vec;
			phi_cam = LzVec_p_rec[num_comb].Vect().Phi();
			
			if(phi_cam > 2*TMath::Pi()){ phi_cam = phi_cam - 2*TMath::Pi();}
                        if(phi_cam < 0){phi_cam = phi_cam + 2*TMath::Pi();}
			pt_cam = LzVec_p_rec[num_comb].Perp();

			double temp_MM = (LzVecPspec-LzVec_p_rec[num_comb]).M2();
			double Z_interp=cam_inst->GetZA(R_vertex,pos_in_ring_a[num_comb],pos_in_ring_b[num_comb],LzVec_p_rec[num_comb].Phi());
			double temp_delta_Z = Z_interp-pos_in_ring_a[num_comb].Z();
			phi_MM = (LzVecPspec).Vect().Phi();
			if(phi_MM > 2*TMath::Pi()){ phi_MM = phi_MM - 2*TMath::Pi();}
			if(phi_MM < 0){phi_MM = phi_MM + 2*TMath::Pi();}
                        pt_MM = (LzVecPspec).Perp();
			
			double temp_dPhi = phi_MM - phi_cam;
			double temp_dpt = pt_MM - pt_cam;

			bool delta_Z_flag = false;
			bool MM_flag = false;
			bool dpt_flag = false;
			bool dphi_flag = false;
			if(std::fabs(temp_delta_Z) < 16) delta_Z_flag = true;
			if(std::fabs(temp_MM) < 0.3) MM_flag = true;
			if(std::fabs(temp_dpt) < 0.3) dpt_flag = true;
			if(std::fabs(temp_dPhi) < 0.4) dphi_flag = true;


			tot_flag = delta_Z_flag && MM_flag && dpt_flag && dphi_flag;

			if(tot_flag){
				//MM_camera = temp_MM;	
				//delta_Z = temp_delta_Z;
				//dPhi = temp_dPhi;
				//dpt = temp_dpt;
				LzVecPro = protons[i].p4;
				count_pro++;
			}


			// Perform kinematic fit of current combination and safe reults (fit results, fitted kinematic ....)
			static Fitter* FitInterface = &(Fitter::GetInstance());  
			
			FitInterface->Init(R_vertex,beam_track,outMu_track, LzVec_p_rec[num_comb], pos_in_ring_a[num_comb], pos_in_ring_b[num_comb]); 
			
			FitInterface->Add_Particle(Trackp1,0.13957);
			FitInterface->Add_Particle(Trackp2,0.13957);
			
			FitInterface->SetupFit();
			
			FitInterface->DoFit();


			TLorentzVector muon_in = *(FitInterface->GetMuonIn()->getCurr4Vec());
			TLorentzVector muon_out = *(FitInterface->GetMuonOut()->getCurr4Vec());
			TLorentzVector proton_target = *(FitInterface->GetProtonTarget()->getCurr4Vec());
			TLorentzVector proton_out = *(FitInterface->GetProtonOut()->getCurr4Vec());
			TLorentzVector pion_1 = *(FitInterface->GetOutParticles()[0]->getCurr4Vec());
			TLorentzVector pion_2 = *(FitInterface->GetOutParticles()[1]->getCurr4Vec());

			


			double temp_chi2;
			int temp_ndf;

			bool fit_converged = FitInterface->GetFitOutput(temp_chi2, temp_ndf);
			
			if(tot_flag && fit_converged){
				LzVecPro = proton_out; //protons[i].p4;
				count_pro++;
			}	
			num_comb++;



		}
			
		if(count_pro != 1) continue;

*/
// reset of variables for next loop
	for (Int_t n=0;n<4 ;n++) {
        	mu_in[n] = 0;
        	mu_out[n] =0;
        	pi_1[n]=0;
        	pi_2[n]= 0;
	//	pro[n] = 0;
  	}

  	for (Int_t n=0;n<4 ;n++) {
       		mu_in[n] = LzVecMu0(n);
      	 	mu_out[n] = LzVecMu1(n);
       		pi_1[n]= LzVecOut1(n);
       		pi_2[n]= LzVecOut2(n);
	//	pro[n] = LzVecPro(n);
  	}

	if (Emiss <= -5 || Emiss >= 20) continue;
	eventFlags.setFlagByName("Emiss_flag", true);
	
	
	
	/*targ_flag = true;
        if(!(target.InTarget(&tempVertex, 1.9))) {
                targ_flag = false;
        }
        if(!target.InTarget(&tempTrack,-318.5,-78.5,1.9)){
                targ_flag = false;
        }
	if(!(targ_flag)){continue;}*/

	if(!(PaAlgo::InTarget(Par_beam,'O',Run,1.9 , 1.2, -318.5, -78.5, 2))){ continue;}
        eventFlags.setFlagByName("vtxInTarget_flag", true);
      	
       	if(!(PaAlgo::CrossCells(beam_track.vTPar(0), Run, 1.9,1.2,-318.5,-78.5,2))){continue;}
	eventFlags.setFlagByName("crossCells_flag", true);	
	
	//if(E_scat > 180 || E_scat < 140){continue;}	
	
	if(E_scat > E){continue;}
	eventFlags.setFlagByName("EE_flag", true);
	if(E > 180 || E < 140){continue;}
	if (!( ((e.TrigMask() & 0x2) != 0) || ((e.TrigMask() & 0x4) != 0) || ((e.TrigMask() & 0x8) != 0) || ((e.TrigMask() & 0x200) != 0) )){
	       //trig_flag = true;
	      continue;
	}
	
	
	eventFlags.setFlagByName("trigger_flag", true);

	
	
	if ((W < 5.0 || W > 17)) continue;
	eventFlags.setFlagByName("W_flag", true);
	if ((Q2 > 10.0) || (Q2 < 1.0)) continue;
	eventFlags.setFlagByName("Q2_flag", true);
	if ((pt2 > 0.5) || (pt2 < 0.01)) continue;
	eventFlags.setFlagByName("pt_flag", true);
	//if ((nu < 20	)) continue;
	//eventFlags.setFlagByName("nu_flag", true);
	if (mom_mv < 15.0) continue; // mom_rho0 cut
	eventFlags.setFlagByName("rhoP_flag", true);
	if ((y > 0.9) || (y < 0.1)) continue;
	eventFlags.setFlagByName("y_flag", true);
	if (M_rho0 < 0.5 || M_rho0 > 1.1) continue; // M_rho0 cut
	eventFlags.setFlagByName("Minv_flag", true);
	if(!TiS_flag) { continue;}
	else{eventFlags.setFlagByName("timeInSpill_flag", true);}

	TVector3 yl1=LzVecOut1.Vect();
  	TVector3 yl2=LzVecOut2.Vect();
  	OpAngle=1000.0*yl1.Angle(yl2);
	tree->Fill();

	e.TagToSave();  

}

 eventFlags.incrementCounters();

}


void UserJobEnd610() {	
 eventFlags.printFlags();
}





