#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "Phast.h"
#include "PaSetup.h"
#include "PaAlgo.h"
#include "PaHodoHelper.h"
#include "PaEvent.h"
#include "G3part.h"

// ******************************************************************************** //
// UserEvent for preselecting exclusive photon events	(DVCS), rho and omega	        //
// In the selection all possible combinations of:					                          //
// DVCS: Vertices (incoming and outgoing muon), exclusive photon and recoil proton	//
// rho: Vertices (incoming and outgoing muon) ... 	//
// omega: Vertices (incoming and outgoing muon) ... 	//
// ******************************************************************************** //

extern "C" float prob_(float&, int&);

static int analysis_flag; 
void userFunction() {
    std::string userInput;

    std::cout << "Please enter the number for your analysis:\n";
    std::cout << "0. Test (DIS) \n";
    std::cout << "1. DVCS\n";
    std::cout << "2. rho\n";
    std::cout << "3. omega\n";
    
    std::cin >> userInput;
    int choice;
    try {
        choice = std::stoi(userInput); // Convert input to integer
    } catch (...) {
        std::cout << "Invalid input. Please enter a number from 0-3.\n";
        return;
    }

    if (choice < 0 || choice > 3) {
        std::cout << "Invalid input. Please enter a number from 0-3.\n";
        return;
    }

    analysis_flag = choice;
    std::cout << "Analysis set to: " << analysis_flag << std::endl;
}

// Global flag for verbose mode 
bool verbose213 = false; // Set to true for verbose output, false to suppress
void printDebug213(const std::string &message) {
  if (verbose213) {
      std::cout << message << std::endl;
  }
};

//*****************************************************************************
void UserEvent213151414(PaEvent & e) {

  // Define constants
  static PaHodoHelper* HodoHelper;

  const double M_pi    = G3partMass[8]; // Pion mass 
  const double M_gamma = G3partMass[1]; // Photon mass
  const double M_mu    = G3partMass[5]; // Muon mass
  const double M_p     = G3partMass[14]; // Proton mass 

  // Declare all objects globally
  static TH1F *h213_cutStats = NULL;
  static TH1F* h213_y     = NULL;
  static TH1F* h213_Q2    = NULL;
  static TH1F* h213_W2    = NULL;
  static TH1F* h213_xbj   = NULL;
  static TH2F* h213_Q2xbj = NULL;
  static TH1F* h213_t     = NULL;
  static TH1F* h213_nu    = NULL;

  static TTree* tree(NULL);

  //
  // Variables to be stored into analysis Ntuple
  //
  // (if you want to extend the list:
  // - add variable here,
  // - add in Ntuple definition below 
  // and do not forget to fill the variable in the code
  //
  static int    Run;         // run number
  static int    Year;        // year 
  static int    Evt;         // event number - unique evt number (based on run, spill and evt num in spill) 
  static int    EvtInSpill;  // event number in spill 
  static int    Spill;       // spill number
  static double TimeInSpill; // time in spill 
  static float  Zprim;       // Z coordinate of primary vertex (10000 if not found)
  static float  Yprim;       // Y coordinate of primary vertex (10000 if not found)
  static float  Xprim;       // X coordinate of primary vertex (10000 if not found)
  static int    trig_mask;

  static double inMu_pz; // Z component of the muon beam momentum 
  static double inMu_py; // Y component of the muon beam momentum 
  static double inMu_px; // X component of the muon beam momentum
  static double inMu_p;  // Magnitude of the muon beam momentum 
  static double inMu_E;  // Magnitude of the muon beam energy

  
  static int    Hodo_i_omu;
  static double outMu_pz; // Z component of the scattered muon momentum 
  static double outMu_py; // Y component of the scattered muon momentum
  static double outMu_px; // X component of the scattered muon momentum 
  static double outMu_p;  // Magnitude of the scattered muon momentum 
  static double outMu_E;  // Magnitude of the scattered muon energy

  static double y;   // fractional energy loss of the incoming lepton 
  static double nu;  // energy of the virtual photon 
  static double Q2;  // four-momentum transfer sqaured (Q is the four-momentum transferred between the incoming and outgoing lepton)
  static double W2;  // effective mass of the final state hadron system squared 
  static double xbj; // measure of the elasticity of the scattering process 
  static double t;   // TODO: add definition 

  // TODO: add flags for other particles 
  // Particle selection flags
  struct OUTMuFlags {
    bool hodo_flag = false;
    bool outMuTrack_flag = false;
    bool vertex_flag = false;
    bool charge_flag = false;
    bool zFirstLast_flag = false;
    bool Q2_flag = false;
    bool y_flag = false;
  };

  struct BeamFlags {
    bool allVtx_flag = false;
    bool pVtx_flag = false;
    bool inMu_flag = false;
    bool inMuTrack_flag = false;
    bool inMuPar_flag = false;
    bool passTarget_flag = false;
    bool zFirst_flag = false;
    bool BMS_flag = false;
    bool FI_flag = false;
    bool SI_flag = false;
    bool momRange_flag = false;
    bool momErr_flag = false;
    bool meantime_flag = false;
    bool TiS_flag = false;
    bool flux_flag = false;
  };

  struct PhotonFlags {
  };

  struct PionFlags {
  };

  struct ProtonFLags {
  };

  struct SpecialFlags {
    bool trig_flag = false;           // ALL 
    bool singleTrack_flag = false;    // DVCS - exclusive cut 
    bool saveDVCSEvent_flag = false;  // DVCS 
    bool saveRhoEvent_flag = false;   // rho 
    bool saveOmegaEvent_flag = false; // omega
  };

  BeamFlags    inMuFlags;    // Create an instance of scattered muon flags
  OUTMuFlags   outMuFlags;   // Create an instance of incoming muon flags
  PhotonFlags  photonFlags;  // Create an instance of exclusive photon flags
  PionFlags    pionFlags;    // Create an instance of pi+/pi- flags 
  ProtonFLags  protonFlags;  // Create an instance of proton flags
  SpecialFlags specialFlags; // Create an instance of analysis related flags 

  // Selection statics (how many events remain after each cut)
  // TODO: NEED TO DECIDE WHAT CUTS WE ARE TRACK 
  static std::map<std::string, std::pair<int, bool>> cutCounts = {
    {"Cut 00", {0, false}}, 
  };

  //*****************************************************************************

  static bool first(true);
  if (first) { // histograms and Ntuples booking block
    Phast::Ref().HistFileDir("UserEvent213151414");
    HodoHelper = &PaHodoHelper::Init("",true);

    userFunction(); 

    // 1D and 2D
    int nCuts = cutCounts.size();
    h213_cutStats = new TH1F("h97_cutStats", "Event Cuts Breakdown; Cut; Number of Events Removed", nCuts, 0, nCuts);

    h213_y     = new TH1F("h213_y", "Fractional Energy Loss of Incoming Muon (y); y; Events", 100, 0, 1);
    h213_Q2    = new TH1F("h213_Q2", "Four-momentum Transfer Squared (Q2); Q2 [GeV2]; Events", 100, 0, 10);
    h213_W2    = new TH1F("h213_W2", "Effective Mass of final state hadrons Squared (W2); W2 [GeV2]; Events; W2 [GeV2]", 100, 0, 350);
    h213_Q2xbj = new TH2F("h213_Q2xbj", "Kinematic Coverage of Dataset; x_bj; Q2 [GeV2]", 100, 0, 1, 100, 0, 100);
    // TODO: add histogram for t and nu 

    const int nBins = 100;       // Number of bins
    double xMin = 1e-3;          // Minimum x value (avoid 0 because log(0) is undefined)
    double xMax = 2.0;           // Maximum x value
    double binEdges[nBins + 1];  // Bin edges array
    for (int i = 0; i <= nBins; ++i) {
      binEdges[i] = xMin * pow(xMax / xMin, double(i) / nBins);
    }
    h213_xbj = new TH1F("h213_xbj", "Elasticity of the Scattering Process (x_bj); x_bj; Events", nBins, binEdges);

    //
    // Ntuple definition 
    //
    tree = new TTree("USR213151414","User 213151414 NTuple"); // name (has to be unique) and title of the Ntuple
    
    tree->Branch("Run",   &Run,   "Run/I");
    tree->Branch("Evt",   &Evt,   "Evt/I");
    tree->Branch("Year",  &Year,  "Year/I");
    tree->Branch("Spill", &Spill, "Spill/I");

    tree->Branch("Zprim",   &Zprim,   "Zprim/F");
    tree->Branch("Yprim",   &Yprim,   "Yprim/F");
    tree->Branch("Xprim",   &Xprim,   "Xprim/F");

    tree->Branch("inMu_pz", &inMu_pz, "inMu_pz/D");
    tree->Branch("inMu_py", &inMu_py, "inMu_py/D");
    tree->Branch("inMu_px", &inMu_px, "inMu_px/D");
    tree->Branch("inMu_E",  &inMu_E,  "inMu_E/D");

    tree->Branch("outMu_pz", &outMu_pz, "outMu_pz/D");
    tree->Branch("outMu_py", &outMu_py, "outMu_py/D");
    tree->Branch("outMu_px", &outMu_px, "outMu_px/D");
    tree->Branch("outMu_E",  &outMu_E,  "outMu_E/D");
    
    // TODO: Add proton (camera), pi+/pi- for rho, exlucsive photon for DVCS, omega ???

    tree->Branch("y",   &y,   "y/D");
    tree->Branch("Q2",  &Q2,  "Q2/D");
    tree->Branch("W2",  &W2,  "W2/D");
    tree->Branch("xbj", &xbj, "xbj/D");
    tree->Branch("t",   &t,   "t/D");
    tree->Branch("nu",  &nu,  "nu/D");

    first=false;
  } // end of histogram booking
  
  //*****************************************************************************
  // Assign names to trigger bits
  enum trigger {   
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
    specialFlags.trig_flag = true;
  }
  if (trig_mask & LT) {
    trigCheck += "LT ";
    specialFlags.trig_flag = true;
  }
  if (trig_mask & OT) {
    trigCheck += "OT ";
    specialFlags.trig_flag = true;
  }
  if (trig_mask & LAST) {
    trigCheck += "LAST ";
    specialFlags.trig_flag = true;
  }

  //*******************************************
  Run        = e.RunNum();
  Year       = e.Year();
  Evt        = e.UniqueEvNum();
  EvtInSpill = e.EvInSpill(); 
  Spill      = e.SpillNum(); 

  // Time in spill check 
  // TODO: ADD CUT CONDIIONS HERE
  TimeInSpill  = e.TimeInSpill(); // check the time in spill 


  // Debug statements ... (1/2)
  printDebug213("      ");
  printDebug213("*** Run: " + std::to_string(Run) + ", spill: " + std::to_string(Spill) + ", event: " + std::to_string(EvtInSpill) + " ***");

  //*******************************************  
  // Reset the "already counted" flags for each event 
  for (auto& cut : cutCounts) {
    cut.second.second = false; 
  }
  cutCounts["Cut 00"].first++;

  // Loop over reconstructed vertices in the event 
  for (int iv = 0; iv < e.NVertex(); iv++) { // begin loop over vertices
  inMuFlags.allVtx_flag = true;

    //******************************************* 
    // Store info about primary vertex (if found) 
    const PaVertex & v = e.vVertex(iv);
    if (!v.IsPrimary()) continue;
    inMuFlags.pVtx_flag = v.IsPrimary() ? true : false;

    //*******************************************
    // Store info about incoming muon beam (inMu)
    int muBeam = v.InParticle(); 
    if (muBeam == -1) continue; // there is incoming particle associated with the primary vertex
    inMuFlags.inMu_flag = (muBeam != -1) ? true : false;  

    const PaParticle & beam = e.vParticle(muBeam);
		int it_beam = beam.iTrack();
		if (it_beam == -1) continue; // the incoming particle has a track associated with it
    inMuFlags.inMuTrack_flag = (it_beam != -1) ? true : false;

    const PaTrack & beam_track = e.vTrack(muBeam);
		if (beam_track.NTPar() == 0) continue; // the track has parameters
    inMuFlags.inMuPar_flag = (beam_track.NTPar() != 0); 

    // PaAlgo::CrossCells(t_beam.vTPar(0),run, Rmax, Ymax, tgt_zmin, tgt_zmax, RmaxMC) - cut 4 starts here
    if (!(PaAlgo::CrossCells(beam_track.vTPar(0),Run, 1.9, 1.2, -318.5, -78.5, 2))) continue;
    inMuFlags.passTarget_flag = PaAlgo::CrossCells(beam_track.vTPar(0), Run, 1.9, 1.2, -318.5, -78.5, 2);

    if (beam_track.ZFirst() >= -78.5) continue; // incoming muon was first measured before the target
    inMuFlags.zFirst_flag = (beam_track.ZFirst() < -78.5) ? true : false;

    // incoming muon is detected by detectors along the beamline  
    int nhits_FI  = beam_track.NHitsFoundInDetect("FI"); 
		int nhits_SI  = beam_track.NHitsFoundInDetect("SI"); 
    int nhits_BMS = beam_track.NHitsFoundInDetect("BM");
    if ((nhits_FI < 2)) continue; 
    inMuFlags.FI_flag = (nhits_FI >= 2);
		if ((nhits_SI < 3)) continue;
    inMuFlags.SI_flag = (nhits_SI >= 3);
    if ((nhits_BMS < 3)) continue; 
    inMuFlags.BMS_flag = (nhits_BMS >= 3); 

    double inMu_mom = beam_track.vTPar(0).Mom();
    if (inMu_mom < 140.0 || inMu_mom > 180.0) continue; // momentum falls within acceptable range
    inMuFlags.momRange_flag = (inMu_mom >= 140.0 && inMu_mom <= 180.0);
    
    double inMu_momErr = sqrt(beam_track.vTPar(0)(5,5))/(beam_track.vTPar(0)(5)*beam_track.vTPar(0)(5));
    if (inMu_momErr > 0.025*inMu_mom) continue; // momentum error falls within acceptable range  
    inMuFlags.momErr_flag = (inMu_momErr <= 0.025*inMu_mom);

    double mean_time = beam_track.MeanTime();  
    if (std::fabs(mean_time) >= 2) continue; // track meantime for incoming muon is within flux requirements 
    inMuFlags.meantime_flag = (std::fabs(mean_time) < 2); 

    // Check if all requirements, beside the time in spill, for proper muon flux are fulfilled
    // TODO: Update to include TiS flag 
    if (inMuFlags.allVtx_flag && inMuFlags.pVtx_flag && inMuFlags.inMu_flag &&
      inMuFlags.inMuTrack_flag && inMuFlags.inMuPar_flag && inMuFlags.passTarget_flag &&
      inMuFlags.zFirst_flag && inMuFlags.FI_flag && inMuFlags.SI_flag 
      && inMuFlags.BMS_flag && inMuFlags.momRange_flag && inMuFlags.momErr_flag &&
      inMuFlags.meantime_flag && !inMuFlags.TiS_flag) {inMuFlags.flux_flag = true;}

    //*******************************************
    // Store info about scattered muon (outMu) 
    // HodoHelper->iMuPrim(v, checkYokeSM2, reject2muEvents, checkCanBeMuon, true, minXX0muPr, true, true) 
    int i_omu = -1; 
		i_omu = HodoHelper->iMuPrim(v,false,false,true,false,15); // index of the scattered muon WITHOUT CHECKING IF IT PASSES the hodoscope check 
		if (i_omu == -1) continue;
    int i_omu_check_hodo = HodoHelper->iMuPrim(v,false,true,true,true,15,true,true); // index for the scattered muon IF IT PASSES the hodoscope check 
    // if scattered muon passed the hodoscope use the corresponding index, if not proceed with other index 
    if (i_omu_check_hodo != -1) {
      outMuFlags.hodo_flag = true;
			i_omu = i_omu_check_hodo;
		}

    const PaParticle & outMu = e.vParticle(i_omu);
    int outMu_itrack = outMu.iTrack();
    if (outMu_itrack == -1) continue; // outgoing muon has a track associated with it 

    const PaTrack & outMu_track = e.vTrack(outMu_itrack); 
    double outMu_mom = outMu_track.vTPar(0).Mom();    

    const PaTPar & Par_beam = beam.ParInVtx(iv); // beam parameters at the vertex
    const PaTPar& Par_outMu = outMu.ParInVtx(iv); // scattered muon parameters at the vertex
    bool isInTarget = PaAlgo::InTarget(Par_beam, 'O', Run, 1.9, 1.2, -318.5, -78.5, 2);
    if (!isInTarget) continue; // vertex is in the target  
    outMuFlags.vertex_flag = isInTarget;

    if (outMu.Q() != beam.Q()) continue; // scattered muon has the same charge as the beam
    outMuFlags.charge_flag = (outMu.Q() == beam.Q());
    
    if (!(outMu_track.ZFirst() < 350 && outMu_track.ZLast() > 350)) continue; // first and last z coordinates are measured before and after SM1
    outMuFlags.zFirstLast_flag = (outMu_track.ZFirst() < 350 && outMu_track.ZLast() > 350);

    //*******************************************
    // Kinematic variables 
    TLorentzVector inMu_TL  = Par_beam.LzVec(M_mu); 
	  TLorentzVector outMu_TL = Par_outMu.LzVec(M_mu); 
    TLorentzVector targ_TL(0,0,0,M_p);
    TLorentzVector q = (inMu_TL - outMu_TL); // four momentum of the virtual photon

    Q2  = PaAlgo::Q2 (inMu_TL, outMu_TL); //Q2  = -(inMu_TL - outMu_TL).M2();
    y   = (inMu_TL.E() - outMu_TL.E()) / inMu_TL.E();
    nu  = (inMu_TL.E() - outMu_TL.E());
    W2  = PaAlgo::W2 (inMu_TL, outMu_TL);
    xbj = PaAlgo::xbj (inMu_TL, outMu_TL); //xbj = Q2/(2*q*targ_TL);

    // TODO: may need to add different kinematic cuts for omega 
    if (analysis_flag != 3) {
        if (Q2 < 1) continue;  // Skip this event if Q2 < 1
        outMuFlags.Q2_flag = (Q2 > 1);

        if (y < 0.01 || y > 0.99) continue; // Skip this event if y is out of range
        outMuFlags.y_flag = (y > 0.01 && y < 0.99);
    }

    //**************************************************
    // ALL CUTS THAT PROCEED DEPEND ON THE ANALYSIS FLAG
    //**************************************************

    //*******************************************
    // Exclusive DVCS selection starts here ... (analysis_flag = 1)
    if (analysis_flag == 1) {
      // Only one outgoing particle (scattered proton and photon are detected using ECals and CAMERA so will not be found here)
      if (v.NOutParticles()!= 1) continue; 
      specialFlags.singleTrack_flag = (v.NOutParticles() == 1); 

      // Check that cut flags are satisfied so far 
      if (inMuFlags.flux_flag && outMuFlags.hodo_flag && outMuFlags.vertex_flag && outMuFlags.charge_flag 
        && outMuFlags.zFirstLast_flag && specialFlags.trig_flag && outMuFlags.Q2_flag && outMuFlags.y_flag) 
        {specialFlags.saveDVCSEvent_flag = true;}
    } 

    //*******************************************
    // Rho selection starts here ... (analysis_flag = 2)
    if (analysis_flag == 2) {
    } 
    
    //*******************************************
    // Debug statements ... (2/2)
    printDebug213("    Vertex: (" + std::to_string(Xprim) + ", " + std::to_string(Yprim) + ", " + std::to_string(Zprim) + ")");
    printDebug213("    mu: P: " + std::to_string(inMu_mom) + " GeV/c, Charge: " + std::to_string(beam.Q()));
    printDebug213("    mu': P: " + std::to_string(outMu_mom) + " GeV/c, Charge: " + std::to_string(outMu.Q()));
    printDebug213("    Kinematics: Q2: " + std::to_string(Q2) +  " GeV2, y: " + std::to_string(y) + ", W2: " + std::to_string(W2) + " GeV2, x: " + std::to_string(xbj));
    printDebug213("    Active triggers: " + trigCheck);
    printDebug213("    Save Event: DVCS" + std::to_string(specialFlags.saveDVCSEvent_flag) +
              ", rho: " + std::to_string(specialFlags.saveRhoEvent_flag) + 
              ", omega: " + std::to_string(specialFlags.saveOmegaEvent_flag)); 

    //*******************************************
    // Fill the tree/histograms with the extracted event information 
    int bin = 1;
    for (const auto &cut : cutCounts) {
        h213_cutStats->SetBinContent(bin, cut.second.first);  // Use cut.second.first to get the count
        h213_cutStats->GetXaxis()->SetBinLabel(bin, cut.first.c_str());  // Set cut labels
        bin++;
    }

    inMu_pz = beam_track.vTPar(0).Pz();
    inMu_py = beam_track.vTPar(0).Py();
    inMu_px = beam_track.vTPar(0).Px();
    inMu_p  = beam_track.vTPar(0).Mom();
    // TODO: get value inMu_E 

    outMu_pz = outMu_track.vTPar(0).Pz();
    outMu_py = outMu_track.vTPar(0).Py();
    outMu_px = outMu_track.vTPar(0).Px();
    outMu_p  = outMu_track.vTPar(0).Mom();
    // TODO get value outMu_E

    tree->Fill();

    // TODO: fill t and nu histograms here as well 
    h213_y->Fill(y);
    h213_Q2->Fill(Q2);
    h213_W2->Fill(W2);
    h213_xbj->Fill(xbj); 
    h213_Q2xbj->Fill(xbj,Q2);

  } // end of loop over vertices

  if (analysis_flag == 0) {
      e.TagToSave();
  }

  if (analysis_flag == 1 && specialFlags.saveDVCSEvent_flag) {
      e.TagToSave();
  }

  if (analysis_flag == 2 && specialFlags.saveRhoEvent_flag) {
      e.TagToSave();
  }

  if (analysis_flag == 3 && specialFlags.saveOmegaEvent_flag) {
    e.TagToSave();
  }

}
