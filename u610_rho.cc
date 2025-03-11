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
#include "TiS_range.cc"
#include "UConn_Tools.h"


bool verbose_mode = true; // Create an instance of verbose_mode

PaHodoHelper* HodoHelper; 
void UserEvent610(PaEvent& e)
{
    //  NLUDATA ld;                   // create NLUDATA structure
    //  vector<LUJET> lujets;         // create LUJET structure
    //  int nparticles;
    //  e.MCgen(ld);
    //  e.MCgen(nparticles,lujets);

  Bool_t intarget;
  static PaHodoHelper* HodoHelper;
 // static TiSRange*     tis_range;
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

  Double_t mu_in[4];
  Double_t mu_out[4];
  Double_t pi_1[4];
  Double_t pi_2[4];

  TLorentzVector LzVecMu0;   // beam mu Lorentz vector
  TLorentzVector LzVecMu1;   // scattered mu Lorentz vector
  TLorentzVector LzVecOut1;  // pi+
  TLorentzVector LzVecOut2;  // pi-
  TLorentzVector LzVecOut;   // pi+ + pi-
//  TLorentzVector LzVecPi0;   // pi0
  TLorentzVector LzVecPspec;

 // int Nout_gen;i
   static TargetCell& target = TargetCell::Instance();
   bool TiS_flag  = false;
   struct EventFlags {
        //DIS cuts
        bool flux_flag  = false;
        bool outMu_flag = false;
        bool Q2_flag    = false;
        bool y_flag     = false;
        bool DIS_flag   = false; // true if all DIS cuts have been applied
        //Exclusive cuts
        bool singleTrack_flag = false;
        bool checkCls_flag    = false; // true if all found clusters have a track, and pass the timing and ECal energy cuts
        bool singleCl_flag    = false; // true after a single neutral cluster has been found in the ECals
        bool exclEvent_flag   = false; // true if all exclusive cuts have ben applied
        bool saveEvent_flag   = false;
    };

    EventFlags eventFlags;   // Create an instance of event flags
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


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
time_in_spill= e.TimeInSpill();
if (Run != LastRun) { // Reinitialize HodoHelper and tis_range only if the run number changes 
       // HodoHelper = & PaHodoHelper::Init("", true);  
        //tis_range  = new TiSRange("/afs/cern.ch/user/n/ntrotta/workDir/work_phast_DVCS2016/phast/user/flux_files");
	Set_TiSrange("/afs/cern.ch/user/n/ntrotta/workDir/work_phast_DVCS2016/phast/user/flux_files/flux_Johannes/2016/flux_files",Run,Run);
	LastRun = Run;  // Update LastRun to the current run number
}
TiS_flag = true; //Check_TiS_window(Run,Spill,time_in_spill);

for (int iv = 0; iv < e.NVertex(); iv++) { // begin loop over vertices
		
	if (!( ((e.TrigMask() & 0x2) != 0) || ((e.TrigMask() & 0x4) != 0) || ((e.TrigMask() & 0x8) != 0) || ((e.TrigMask() & 0x200) != 0) )) continue;
	const PaVertex & v = e.vVertex(iv);
	Nout = v.NOutParticles(); // number of tracks in vertex
	
	if(Nout !=3) continue;
	if (!v.IsPrimary()) continue;


	Zprim = v.Pos(2);
	Yprim = v.Pos(1);
	Xprim = v.Pos(0);
	Nout = v.NOutParticles(); // number of tracks in vertex
	
	if(Nout !=3) continue;
	static PaParticle beam; 
      	static PaTrack beam_track; 
      	static PaTPar Par_beam;

	static BeamFluxParams beamParams; // Create an instance of BeamFluxParams
	eventFlags.flux_flag = beamFluxCheck(e, v, iv, Run, TiS_flag, beamParams, beam, beam_track, Par_beam);
	if (!eventFlags.flux_flag) continue;
	
	static PaParticle outMu; 
      	static PaTrack outMu_track; 
      	static PaTPar Par_outMu;

      	static OutMuParams outMuParams; // Create an instance of OutMuParams 
      	eventFlags.outMu_flag = outMuCheck(e, v, iv, Run, beam, HodoHelper, true, outMuParams, outMu, outMu_track, Par_outMu);
      	if (!eventFlags.outMu_flag) continue;
	
	int imu1 = v.iMuPrim(false,true,true,false,15);	
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
       
    if(Trackp1.XX0() > 10) continue; // penetration length of hadron track shorter than 10 radiation lengths
    
    if(Trackp2.XX0() > 10) continue; // penetration length of hadron track shorter than 10 radiation lengths

    if(Trackp1.Chi2tot()/Trackp1.Ndf() > 10) continue; // good fit quality of outgoing hadron ch_red^2 < 10
       
    if(Trackp2.Chi2tot()/Trackp2.Ndf() > 10) continue; // good fit quality of outgoing hadron ch_red^2 < 10
       
    if(Trackp1.ZFirst()>350 || Trackp2.ZFirst()>350) continue;// tracks of hadrons start upstream from SM1
                 
    if(h1.Q() == h2.Q()) continue; //skip equaly charged particles 

       if(h1.Q() < h2.Q()) swap(h1,h2); //now first particle is positive            
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

// reset of variables for next loop
	for (Int_t n=0;n<4 ;n++) {
        	mu_in[n] = 0;
        	mu_out[n] =0;
        	pi_1[n]=0;
        	pi_2[n]= 0;
  	}

  	for (Int_t n=0;n<4 ;n++) {
       		mu_in[n] = LzVecMu0(n);
      	 	mu_out[n] = LzVecMu1(n);
       		pi_1[n]= LzVecOut1(n);
       		pi_2[n]= LzVecOut2(n);
  	}

	TVector3 yl1=LzVecOut1.Vect();
  	TVector3 yl2=LzVecOut2.Vect();
  	OpAngle=1000.0*yl1.Angle(yl2);


	tree->Fill();

	e.TagToSave();  

} 
	





}
