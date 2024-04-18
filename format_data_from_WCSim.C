/**
 * @file get_sig_from_WCSim.C
 * @author Antoine Beauchêne, cleaned and upgraded by Lorenzo Perisse
 * @brief This macro creates a tree from simulated signal (ie with no noise) *
 * @version 1.1.0
 * @date May 2023
 * @copyright Copyright (c) 2022
 *
 */



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <iterator>
#include <vector>
#include <cstdlib>
#include <dirent.h>
#include <unistd.h>
#include <unordered_map>
#include <numeric>

using namespace std;

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TRandom1.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TBox.h>
#include <TLine.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TMath.h>
#include <TObject.h>
#include <TGraph.h>

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimPmtInfo.hh"


#define VERBOSE 
#define PRINT_EVERY 1000        //< Print info every PRINT_EVERY events if VERBOSE is defined

#define MAX_HITS 1000000          //< Max number of hits recorded during 1ms of noise with 4.2kHz of darknoise (depends on the darknoise frequency)

#define TIME_GAP 2000.            //< Time gap between event insertion
#define TIME_DURATION_BG 1e6      //< Time of noise events
#define TIME_WINDOW 500.          //< Time of the window used in the SWA. HK=400ns, SK=200ns,

#define DRAW_TIME 20000.          //< Upper time limit in ns for the plots
#define DRAW_HITY 30.             //< Upper number of hits



struct HitInfo {
  float X;
  float Y;
  float Z;
  float T;
  float Charge;
  int Id;
  int Type;
  int EventNumber;
};


struct TimeWindow {
    int Index;
    int Hits;
    int Type;
    float TriggerTime;
};


struct arguments
{
  string inputNoise  = "";
  string inputSignal = "";
  bool signalRate = false;
  int eventType = -1;
  int PMT = 19547;             //HK value
  float tankRadius = 3240.0;   //HK value
  float tankHalfz = 3287.55;   //HK value
  float darkRate = 4.2;
  int windowSignal = 1;
  int nhits200ns = 25;
  int maxEventBg = -1;
  bool verbose = false;
};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

arguments FetchInput(int argc, char* argv[])
{
    int c = -1;
    arguments arglist;

    //Input in c the argument (-f etc...) and in optarg the next argument
    //When the above test becomes -1, it means it fails to find a new argument
    cout << endl;
    while( (c = getopt(argc, argv, "n:s:e:r:p:t:z:d:w:h:m:v:o")) != -1 )
    {
        switch(c)
        {
        //Input file name for noise
        case 'n':
            arglist.inputNoise = optarg;
            cout << "Input WCSim file for noise: " << arglist.inputNoise << endl;
            break; 

        //Input file name for signal
        case 's':
            arglist.inputSignal = optarg;
            cout << "Input WCSim file for signal: " << arglist.inputSignal << endl;
            break; 

        //Input signal event rate
        case 'r':
            arglist.signalRate = optarg;
            if(arglist.signalRate){cout << "Input signal event rate: " << arglist.signalRate << endl;}
            else                  {cout << "Input signal event rate: 1 event every " << TIME_GAP << " ns" << endl;}
            break;

        //Event type
        case 'e':
            arglist.eventType = atoi(optarg);
            cout << "Event labeled as type #" << arglist.eventType << endl;
            break;

        //Number of PMT in the considered detector
        case 'p':
            arglist.PMT = atoi(optarg);
            cout << "Number of PMT in the considered detector: " << arglist.PMT << endl;
            break;

       //Tank radius
        case 't':
            arglist.tankRadius = atof(optarg);
            cout << "Tank radius: " << arglist.tankRadius << " cm" << endl;
            break;

       //Tank half height
        case 'z':
            arglist.tankHalfz = atof(optarg);
            cout << "Tank half height: " << arglist.tankHalfz << " cm" << endl;
            break;

        //Darkrate in kHz
        case 'd':
            arglist.darkRate = atof(optarg);
            cout << "Dark rate frequency: " << arglist.darkRate << " kHz" << endl;
            break;

        //Minimum hit of signal per trigger window to define the window as signal
        case 'w':
            arglist.windowSignal = atoi(optarg);
            cout << "Minimum hit of signal per trigger window to define the window as signal: " << arglist.windowSignal << endl;
            break;

        //Number of hits during 200ns needed to issue a trigger in addition to the average noise rate
        case 'h':
            arglist.nhits200ns = atoi(optarg);
            cout << "Number of hits during 200ns needed to issue a trigger in addition to the average noise rate: " << arglist.nhits200ns << endl;
            break;

        //Number of signal event to treat
        case 'm':
            arglist.maxEventBg = atoi(optarg);
            cout << "Number of background events to extract: " << arglist.maxEventBg << endl;
            break;

        //Warning
        case 'v':
            arglist.verbose = true;
            cout << "VERBOSE option on" << endl;
            break;
        }
    }
    cout << endl;

    return arglist;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Trigger algorithm
// Based on the algorithm used in WCSim
// See WCSimWCTriggerBase::AlgNDigit() from WCSimWCTrigger.cc and WCSimWCTrigger.hh


// Find triggers
vector<TimeWindow> AlgTrigger(int nhits, vector<float> *digit_times, int ndigit_per_200ns=25, int npmts=19547, float darkrate=4.2)
{
  //Sanity check
  if(digit_times->size() != nhits){cout << "ERROR, number of digihits 1= size of digihit time vector..." << endl;}

  float dark_rate_Hz = darkrate * 1000.;                                      //Dark rate of PMTs in Hz
  float trigger_window_seconds = TIME_WINDOW * 1E-9;                          //Length of the trigger window in seconds
  float average_occupancy = dark_rate_Hz * trigger_window_seconds * npmts;    //Average number of dark rate hits in a trigger window
  float size_window  = TIME_WINDOW;                                           //Size of the trigger window in ns
  int ndigitsThreshold = round(average_occupancy  +  ndigit_per_200ns*TIME_WINDOW/200.);    //Number of hits in a trigger window to issue trigger (for instance for SK: ndigitsThreshold = 25 for size_window = 200)
  
  TimeWindow window_tmp;
             window_tmp.Index = 0;
             window_tmp.Hits  = 0;
             window_tmp.Type  = 0;
             window_tmp.TriggerTime = 0.;
  vector<TimeWindow> tabWin;                  //Table of index of the first hit ([0]), number of hits ([1]) in the windows, and if the type of the window ([3]: -1 for noise only, event # otherwise). Each slot can potentially hold a trigger window
  vector<TimeWindow> all_possible_tabWin;     //List of all the possible time window, ie one time window starting at each hit until the end of the file
  vector<TimeWindow> all_local_max_tabWin;    //List of all the possible time window, ie one time window starting at each hit until the end of the file


  //Get the range of digit_times to loop over (equivalent of looping over PMTs and Digithits in each PMT).  If ndigits > Threshhold in a trigger window, then we have a trigger
  int ntrig = 0;
  float firsthit = (*digit_times)[0];
  float lasthit  = (*digit_times)[nhits-1];
  float window_step_size = 5; //step the search window along this amount if no trigger is found
  float window_start_time  = firsthit;
        window_start_time -= int(window_start_time) % 5;
  float window_end_time    = lasthit - size_window + window_step_size;
        window_end_time   += int(window_end_time) % 5;
  float window_start_time_saveproof = window_start_time;


  // Store all the possible time windows ; the upper time limit is set to the final possible full trigger window
  while(window_start_time <= window_end_time)
  {
    window_tmp.Index = 0;
    window_tmp.Hits  = 0;
    window_tmp.Type  = 0;
    window_tmp.TriggerTime = 0.;

    int n_digits = 0;
    int first_hit_in_window = 0;
    float next_hit_time = window_end_time;
    float triggertime; //trigger time is the time of the first hit above threshold
    bool IsFirstHitInWindow = false;
    
    //Loop over sorted digithits until you issue a trigger
    for (int i = 0 ; i < nhits ; i++)
    {
      float digit_time = (*digit_times)[i];

      //Is the i^{th} hit the first hit of the trigger window
      if (digit_time>=window_start_time  &&  IsFirstHitInWindow==false)
      {
        first_hit_in_window = i;
        next_hit_time = (*digit_times)[i+1];    //Time of the next hit
        IsFirstHitInWindow = true;
      }

      //Is the i^{th} hit in trigger window?
      if (digit_time >= window_start_time   &&   digit_time <= (window_start_time + size_window)){n_digits++;}
    }

    window_tmp.Index = first_hit_in_window;
    window_tmp.Hits  = n_digits;
    window_tmp.Type  = 0;
    window_tmp.TriggerTime = 0.;

    // Store the possible time windows
    all_possible_tabWin.push_back(window_tmp);

    //Advance by one time step (5 ns) ; it must be at least bring you to the next hit, otherwise the window will be the same as at the previous step
    window_start_time += window_step_size;
    if (window_start_time < next_hit_time)
    {
      window_start_time  = next_hit_time;
      window_start_time -= int(window_start_time) % 5;
    }

    //Force out of the loop when...
    if(window_start_time + size_window >= window_end_time){break;}      //... windows goes out of the timeframe
    if  (window_start_time_saveproof > window_start_time){break;}         //... windows goes back in time
    else{window_start_time_saveproof = window_start_time;}
  }



  // Now, we loop over all the possible time windows and we keep those with a number of hits > threshold
  for (int i=all_possible_tabWin.size()-1 ; i>=0 ; --i)
  {
      if (all_possible_tabWin.at(i).Hits < ndigitsThreshold) 
      {
          all_possible_tabWin.erase(all_possible_tabWin.begin() + i);
      }
  }


  // Now, we loop over all the local maximum time windows and we remove those that overlap and have lower number of hits ; then we issue trigger
  int index_check = 0;
  for(int i=0 ; i<all_possible_tabWin.size()-1 ; i++)
  {
    window_tmp.Index = 0;
    window_tmp.Hits  = 0;
    window_tmp.Type  = 0;
    window_tmp.TriggerTime = 0.;

    // The 2 consecutives windows do not overlap, then save the first and then go to next comparison
    if(all_possible_tabWin.at(index_check).Index + all_possible_tabWin.at(index_check).Hits -1  <  all_possible_tabWin.at(i+1).Index)
    {
      tabWin.push_back( all_possible_tabWin.at(index_check) );
      index_check = i+1;
      ntrig++;
    }

    // The 2 consecutives windows overlap
    else if(all_possible_tabWin.at(index_check).Index + all_possible_tabWin.at(index_check).Hits -1  >=  all_possible_tabWin.at(i+1).Index)
    {
      if(all_possible_tabWin.at(index_check).Hits <= all_possible_tabWin.at(i+1).Hits){index_check++;}    // 2nd window has more hits, so keep the 2nd window and proceed to next iteration
    }
  }
  

#ifdef VERBOSE
  cout << "AlgTrigger. Found first/last hits. Looping from " << firsthit << " ns to " << window_end_time << " ns in steps of " << window_step_size << " ns" << endl;
  cout << "AlgTrigger. Number of entries in input digithit collection: " << nhits << endl;
  cout << "Trigger hit threshold: " << ndigitsThreshold << " hits" << endl;
  cout << "Found " << ntrig << " triggers" << endl;
#endif

  return tabWin;
}







////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Get the hits information of events from a WCSim file
 *
 * @param infilename Name of the file containing the WCSim events
 * @param outfilename Name of the ROOT file in which the tree containing the relevant information of each event
 */

void get_hits_from_WCSim(string infilename, string outfilename, int event_type)
{
  int Ipnu = 0;
  int ntrack=0; 
  int nevents=0, nevents_check=0;
  double dr, dx, dy, dz, dt;
  double event_duration = 0.;
  double total_event_duration = 0.;
  bool IsSignal = false;
  TRandom *r1 = new TRandom1();

  // Declare variables used to store data from the branches of the event TTree
  int n_hits;                                           //< Number of hits in event
  int eventType;                                    //< Type of the event
  int hitsType[MAX_HITS];                               //< Indicate if the hit is from noise (0) or from the event (1)
  int tubeIds[MAX_HITS];                                //< tubeID of each PMT registering a photo-electron
  float vertex_x, vertex_y, vertex_z, vertex_t;                                      //< x, y, z and t coordinates of interaction vertex
  float direction_x, direction_y, direction_z; //< x, y, z arrays of hit coordinates
  float energy;                                         //< Energy of the incoming particle
  float hitx[MAX_HITS], hity[MAX_HITS], hitz[MAX_HITS], hitc[MAX_HITS], hitt[MAX_HITS]; //< x, y, z, charge and time  arrays of hit coordinates


/*************************************************************************/
//Prepare data

  //Open the file containing the tree.
  TFile *file = TFile::Open(infilename.data(), "READ");
  TFile *outfile = new TFile(outfilename.data(), "RECREATE");
  TTree *outtree = new TTree("THits", "PMT hits from event");

  //Creates the TTree object, to read detector geometry data, contains detector info - only need 1 "event"
  TTree *geotree=(TTree*)file->Get("wcsimGeoT");
  WCSimRootGeom *geo = new WCSimRootGeom();
  geotree->SetBranchAddress("wcsimrootgeom", &geo);
  if(geotree->GetEntries() == 0){cout << "ERROR, there is 0 event in the Geometry TTree" << endl; exit(2);}
  geotree->GetEntry(0);   //Get event from geometry TTree

  //Creates the TTree object, to read event data
	TTree *wcsimT=(TTree*)file->Get("wcsimT"); 					    
  WCSimRootEvent *wcsimrootsuperevent = new WCSimRootEvent();         // Create a WCSimRootEvent to store data from the tree in

  // Set the branch address for reading from the trees
  TBranch *branch;
  branch = wcsimT->GetBranch("wcsimrootevent");
  branch->SetAddress(&wcsimrootsuperevent);
  wcsimT->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);     // Force deletion to prevent memory leak 
  nevents = wcsimT->GetEntries();

  // Make tree branches
  outtree->Branch("eventType", &eventType, "eventType/I");
  outtree->Branch("vertex_x", &vertex_x, "vertex_x/F");
  outtree->Branch("vertex_y", &vertex_y, "vertex_y/F");
  outtree->Branch("vertex_z", &vertex_z, "vertex_z/F");
  outtree->Branch("vertex_t", &vertex_t, "vertex_t/F");
  outtree->Branch("direction_x", &direction_x, "direction_x/F");
  outtree->Branch("direction_y", &direction_y, "direction_y/F");
  outtree->Branch("direction_z", &direction_z, "direction_z/F");
  outtree->Branch("energy", &energy, "energy/F");
  outtree->Branch("n_hits", &n_hits, "n_hits/I");
  outtree->Branch("hitsType", &hitsType, "hitsType[n_hits]/I");
  outtree->Branch("hitx", hitx, "hitx[n_hits]/F");
  outtree->Branch("hity", hity, "hity[n_hits]/F");
  outtree->Branch("hitz", hitz, "hitz[n_hits]/F");
  outtree->Branch("hitt", hitt, "hitt[n_hits]/F");
  outtree->Branch("hitc", hitc, "hitc[n_hits]/F");
  outtree->Branch("tubeIds", tubeIds, "tubeIds[n_hits]/I");


/*************************************************************************/
// Loop over the events

  eventType = 0;
  for(int i=0 ; i<nevents ; i++)
  {

    wcsimT->GetEvent(i);                                                    // Load the i^th event into wcsimrootsuperevent
    WCSimRootTrigger *wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);    // Load the first (and only) trigger of this event into wcsimrootevent

    // For each event, loop over the associated hits (at least 1!) and fill the branches with hits coordinates, charge and time of each hit
    n_hits = wcsimrootevent->GetNcherenkovdigihits();
    if(n_hits >= 1)
    {

      /***************************************/
      /*                                     */
      /* THIS BLOC IS ONLY FOR SIGNAL EVENTS */
      /*                                     */
      /***************************************/

      //Event vertex
      vertex_x = wcsimrootevent->GetVtx(0);
      vertex_y = wcsimrootevent->GetVtx(1);
      vertex_z = wcsimrootevent->GetVtx(2);

      // For each event, loop over the associated tracks and get back the energy
      ntrack = wcsimrootevent->GetNtrack();
      for(int k=0 ; k<ntrack ; k++)
      {
        TObject *element = (wcsimrootevent->GetTracks())->At(k);
        WCSimRootTrack *wcsimroottrack = (WCSimRootTrack*) (element);

        // Ipnu: see PDG particle numbering scheme (https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf or https://particle.wiki/wiki/PDG_particle_numbering_scheme)
        //  Flag: -1 = incoming particle ; -2 = target ; 1 = outgoing lepton ; 2 = most energetic outgoing nucleon
        if(wcsimroottrack->GetIpnu()!=0  &&  wcsimroottrack->GetFlag()==-1)
        {
          Ipnu = wcsimroottrack->GetIpnu();
          eventType   = event_type;
          energy      = wcsimroottrack->GetE();
          vertex_t    = wcsimroottrack->GetTime();
          direction_x = wcsimroottrack->GetDir(0);
          direction_y = wcsimroottrack->GetDir(1);
          direction_z = wcsimroottrack->GetDir(2);
          IsSignal    = true;      // This is a signal file
        }
      }



      /*************************************************/
      /*                                               */
      /* THIS BLOC IS FOR BOTH SIGNAL AND NOISE EVENTS */
      /*                                               */
      /*************************************************/

      for(int j=0 ; j<n_hits ; j++)
      {
        TObject *Hit; Hit = (wcsimrootevent->GetCherenkovDigiHits())->At(j);
        
        WCSimRootCherenkovDigiHit *cDigiHit = dynamic_cast<WCSimRootCherenkovDigiHit*>(Hit);
        if(IsSignal == true){hitsType[j]= 1;}
        else                {hitsType[j]= 0;}
        tubeIds[j] = cDigiHit->GetTubeId();
        hitt[j]    = cDigiHit->GetT();
        hitc[j]    = cDigiHit->GetQ();

        // WCSimRootPMT hitpmt; hitpmt = geo->GetPMT(tubeIds[j]-1, false);     //For WCSim-1.12 and above
        WCSimRootPMT hitpmt; hitpmt = geo->GetPMT(tubeIds[j]-1);     //For WCSim-1.10 and below
        hitx[j]  = hitpmt.GetPosition(0);
        hity[j]  = hitpmt.GetPosition(1);
        hitz[j]  = hitpmt.GetPosition(2);
      }


      outtree->Fill();    //Fill only if there is at least one hit
      nevents_check++;
      event_duration = *max_element(hitt , hitt + n_hits);
      total_event_duration += event_duration;

  #ifdef VERBOSE
      if (i % PRINT_EVERY == 0)
      {
        cout << "Event #" << i << " / " << nevents
             << ",  Cherenkov hits=" << n_hits;
        if (ntrack > 0)
        {
        cout << ",  particle ID=" << Ipnu
             << ",  Total energy=" << energy
             << ",  t0=" << vertex_t 
             << " ns up to t1=" << event_duration << " ns";
        }
        cout << endl;
      }
#endif

    }
  }

  outtree->Write();
  outfile->Close();


/*************************************************************************/
// Print

  cout << "There are " << nevents       << " events in the input file "  << infilename  << endl;
  cout << "There are " << nevents_check << " events in the output file " << outfilename << endl;
  cout << "The total event duration in the output file is " << total_event_duration << " ns.\n\n" << endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief This function merges the signal hits contained in `sigfilename` with the noise windows contained in `noisefilename` and stores the data in `outfilename`
 *
 * @param sigfilename Name of the file containing the simulated event data, the signal
 * @param noisefilename Name of the file containing the 1ms noise windows data
 * @param outfilename Name of the file containing the data from the signal and noise merging
 */

int merge_signal_to_noise(string sigfilename, string noisefilename, string outfilename, double time_gap)
{
  // The suffix "_SIG" is used for data related to SIGnal file
  int n_hits_SIG;
  float energy_SIG;
  float vertex_x_SIG;
  float vertex_y_SIG;
  float vertex_z_SIG;
  float vertex_t_SIG;
  float direction_x_SIG;
  float direction_y_SIG;
  float direction_z_SIG;
  int hitsType_SIG[MAX_HITS];
  int tubeIds_SIG[MAX_HITS];
  float hitx_SIG[MAX_HITS];
  float hity_SIG[MAX_HITS];
  float hitz_SIG[MAX_HITS];
  float hitt_SIG[MAX_HITS];
  float hitc_SIG[MAX_HITS];

  // The suffix "_BG" is used for data related to noise file (BackGround)
  int n_hits_BG;
  int hitsType_BG[MAX_HITS];
  int tubeIds_BG[MAX_HITS];
  float hitx_BG[MAX_HITS];
  float hity_BG[MAX_HITS];
  float hitz_BG[MAX_HITS];
  float hitt_BG[MAX_HITS];
  float hitc_BG[MAX_HITS];

  // Declare variables for the branches of the output TTree
  int n_hits;                      //< Number of hits in event
  vector<int> *eventNumber = 0;    //< Event number used to differentiate between hits
  vector<float> *hitx = 0;         //< Hits x coordinates
  vector<float> *hity = 0;         //< Hits y coordinates
  vector<float> *hitz = 0;         //< Hits z coordinates
  vector<float> *hitt = 0;         //< Hits timings
  vector<float> *hitc = 0;         //< Hits charges
  vector<int> *tubeIds = 0;
  vector<int> *hitsType = 0;
  vector<float> *vertex_x = 0;
  vector<float> *vertex_y = 0;
  vector<float> *vertex_z = 0;
  vector<float> *vertex_t = 0;
  vector<float> *direction_x = 0;
  vector<float> *direction_y = 0;
  vector<float> *direction_z = 0;
  vector<float> *energy = 0;

  //Variables for event insertion
  int neventsNoise=0, neventsSignal=0;
  int sig_for_one_noise = 0;
  int pool_of_signal = 0;
  int first_nb_of_entries=0, second_nb_of_entries=0;
  int index_SIG=0, index_BG=0;
  int pool_SIG=0, pool_BG=0;
  float time_current=0., time_wall=0., time_limit=0.;
  float max_hit_signal_time = 0.;
  TRandom  *r1 = new TRandom1(0);   // Use random for the inserting time
  TRandom3 *r3 = new TRandom3(0);   // Use random for the inserting time


  TFile *file_SIG = TFile::Open(sigfilename.data(), "READ");
  TTree *tree_SIG = (TTree*)file_SIG->Get("THits"); 	
         tree_SIG->SetBranchAddress("vertex_x", &vertex_x_SIG);
         tree_SIG->SetBranchAddress("vertex_y", &vertex_y_SIG);
         tree_SIG->SetBranchAddress("vertex_z", &vertex_z_SIG);
         tree_SIG->SetBranchAddress("vertex_t", &vertex_t_SIG);
         tree_SIG->SetBranchAddress("direction_x", &direction_x_SIG);
         tree_SIG->SetBranchAddress("direction_y", &direction_y_SIG);
         tree_SIG->SetBranchAddress("direction_z", &direction_z_SIG);
         tree_SIG->SetBranchAddress("energy", &energy_SIG);
         tree_SIG->SetBranchAddress("n_hits", &n_hits_SIG);
         tree_SIG->SetBranchAddress("hitsType", &hitsType_SIG);
         tree_SIG->SetBranchAddress("hitx", &hitx_SIG);
         tree_SIG->SetBranchAddress("hity", &hity_SIG);
         tree_SIG->SetBranchAddress("hitz", &hitz_SIG);
         tree_SIG->SetBranchAddress("hitt", &hitt_SIG);
         tree_SIG->SetBranchAddress("hitc", &hitc_SIG);
         tree_SIG->SetBranchAddress("tubeIds", &tubeIds_SIG);


  TFile *file_BG  = TFile::Open(noisefilename.data(), "READ");
  TTree *tree_BG  = (TTree*)file_BG->Get("THits");
         tree_BG->SetBranchAddress("n_hits", &n_hits_BG);
         tree_BG->SetBranchAddress("hitsType", &hitsType_BG);
         tree_BG->SetBranchAddress("hitx", &hitx_BG);
         tree_BG->SetBranchAddress("hity", &hity_BG);
         tree_BG->SetBranchAddress("hitz", &hitz_BG);
         tree_BG->SetBranchAddress("hitt", &hitt_BG);
         tree_BG->SetBranchAddress("hitc", &hitc_BG);
         tree_BG->SetBranchAddress("tubeIds", &tubeIds_BG);


  // Prepare output file and output tree - 1 outfile only with the same data
  TFile *outfile = new TFile(outfilename.data(), "recreate");
  TTree *outtree = new TTree("THits", "PMT hits from event + background");
         //Vertex
         outtree->Branch("eventNumber", &eventNumber);
         outtree->Branch("vertex_x", &vertex_x);
         outtree->Branch("vertex_y", &vertex_y);
         outtree->Branch("vertex_z", &vertex_z);
         outtree->Branch("vertex_t", &vertex_t);
         outtree->Branch("direction_x", &direction_x);
         outtree->Branch("direction_y", &direction_y);
         outtree->Branch("direction_z", &direction_z);
         outtree->Branch("energy", &energy);
        //Hits
         outtree->Branch("n_hits", &n_hits);
         outtree->Branch("hitsType", &hitsType);
         outtree->Branch("hitx", &hitx);
         outtree->Branch("hity", &hity);
         outtree->Branch("hitz", &hitz);
         outtree->Branch("hitt", &hitt);
         outtree->Branch("hitc", &hitc);
         outtree->Branch("tubeIds", &tubeIds);



/*************************************************************************/
//Insert signal events over noise events
//Take a noise file and paste signal over it with a given rate
//When you reach the end of the noise event, take the next noise event
//Ends when you run up of signal or noise events

  neventsNoise  = tree_BG->GetEntries();                     // Number of noise events ~ 1000 ms
  neventsSignal = tree_SIG->GetEntries();                    // Number of signal events ~ 100 ns
  cout << "Inserting " << neventsSignal << " signal events spaced by " << time_gap << " ns on top of " << neventsNoise << " noise events" << endl;

  // First phase: Loop over all the hits of the signal and the noise ; Stops when there is no more signal or no more noise
  while(pool_BG < neventsNoise   &&   pool_SIG < neventsSignal)
  {
    /**************************/
    /*                        */
    /* THIS BLOC IS FOR NOISE */
    /*                        */
    /**************************/

    tree_BG->GetEntry(pool_BG);
    n_hits  = n_hits_BG;              //Total number of hits of noise during this noise event

    //Fill data of each noise hit
    index_BG = 0;
    while (index_BG < n_hits_BG)
    {
      eventNumber->push_back(-1);     // Hits induced by noise are tagged "-1"
      hitx->push_back(hitx_BG[index_BG]);
      hity->push_back(hity_BG[index_BG]);
      hitz->push_back(hitz_BG[index_BG]);
      hitt->push_back(hitt_BG[index_BG]);
      hitc->push_back(hitc_BG[index_BG]);
      tubeIds->push_back(tubeIds_BG[index_BG]);
      hitsType->push_back(hitsType_BG[index_BG]);
      index_BG++;   //We count this hit as done and go to the next
    }

    //Set absolute time origin and limit
    time_current = *min_element(hitt->begin(), hitt->end());    //Set the atarting time at the minimum of the noise event 
    time_wall    = *max_element(hitt->begin(), hitt->end());
    time_limit   = 0.99*time_wall;                              //We stop the merging when we reach the time_limit ; The last 1% of the noise event should be large enough time-wise to fill at one last event
    cout << "Noise event #" << pool_BG << " completed, starting at t=" << time_current << " ns up to t=" << time_wall << " ns" << endl;
    pool_BG++;      //We count this event as done and go to the next



    /***************************/
    /*                         */
    /* THIS BLOC IS FOR SIGNAL */
    /*                         */
    /***************************/

    // Fill the info of signal
    while(pool_SIG < neventsSignal   &&   time_current < time_limit)
    {
      tree_SIG->GetEntry(pool_SIG);
      if (n_hits_SIG > 0)
      {
        n_hits += n_hits_SIG;                        //Increment by the number of hit due the j^th signal
        time_current += time_gap;                    //Time increase in ns to paste signal further ahead in time (fixed)
        if(time_current > time_limit){break;}        //Stop the merging for this noise event if we go over the noise event max time

        //Fill data of each signal hit
        index_SIG = 0;
        max_hit_signal_time = 0.;
        while (index_SIG < n_hits_SIG)
        {
          eventNumber->push_back(pool_SIG);     // Hits induced by signal are tagged by a number "pool_SIG"
          hitx->push_back(hitx_SIG[index_SIG]);
          hity->push_back(hity_SIG[index_SIG]);
          hitz->push_back(hitz_SIG[index_SIG]);
          hitt->push_back(hitt_SIG[index_SIG] + time_current);     //< Time at which we will insert the signal hit, + time_current because signal time should starts at 0 in WCSim
          hitc->push_back(hitc_SIG[index_SIG]);
          tubeIds->push_back(tubeIds_SIG[index_SIG]);
          hitsType->push_back(hitsType_SIG[index_SIG]);
          index_SIG++;
          if(hitt_SIG[index_SIG] > max_hit_signal_time){max_hit_signal_time = hitt_SIG[index_SIG];}
        }

        vertex_x->push_back(vertex_x_SIG);
        vertex_y->push_back(vertex_y_SIG);
        vertex_z->push_back(vertex_z_SIG);
        vertex_t->push_back(vertex_t_SIG + time_current);
        direction_x->push_back(direction_x_SIG);
        direction_y->push_back(direction_y_SIG);
        direction_z->push_back(direction_z_SIG);
        energy->push_back(energy_SIG);

#ifdef VERBOSE
        if (pool_SIG % PRINT_EVERY == 0){cout << "Signal event #" << pool_SIG << " completed, starting at t=" << time_current << " ns up to t=" << max_hit_signal_time + time_current << " ns (lasting " << max_hit_signal_time << " ns), " << n_hits_SIG << " Cherenkov hits, Total energy " << energy_SIG << " MeV" << endl;}
#endif

        time_current = *max_element(vertex_t->begin(), vertex_t->end());         //Set absolute current time
        if(time_current > time_wall){cout << "ERROR, you went over the time allocated to the noise event (see l.594)..." << endl;  exit(9);} //Stop the merging for this noise event if we go over the noise event max time
      }
      pool_SIG++;                                                       //We count this event as done and go to the next
    }
    n_hits = hitt->size();     //Resize this vector to fit exactly the number of hits


    // Sort all the hit info based on the hit time
    vector<HitInfo> HitInfo_vec;
    for(int i=0 ; i<hitt->size() ; i++)
    {
      HitInfo_vec.push_back({(*hitx)[i], 
                             (*hity)[i], 
                             (*hitz)[i], 
                             (*hitt)[i], 
                             (*hitc)[i], 
                             (*tubeIds)[i], 
                             (*hitsType)[i], 
                             (*eventNumber)[i]
                             });
    }
  	sort(HitInfo_vec.begin(), HitInfo_vec.end(), []( HitInfo i,  HitInfo j) { return i.T < j.T; });
    for(int i=0 ; i<hitt->size() ; i++)
    {
      (*hitx)[i]        = HitInfo_vec[i].X;
      (*hity)[i]        = HitInfo_vec[i].Y;
      (*hitz)[i]        = HitInfo_vec[i].Z;
      (*hitt)[i]        = HitInfo_vec[i].T;
      (*hitc)[i]        = HitInfo_vec[i].Charge;
      (*tubeIds)[i]     = HitInfo_vec[i].Id;
      (*hitsType)[i]    = HitInfo_vec[i].Type;
      (*eventNumber)[i] = HitInfo_vec[i].EventNumber;
    }


    //Fill TTree, reinitialize data vectors and go to next noise event
    outtree->Fill();
    eventNumber->clear();
    hitx->clear();
    hity->clear();
    hitz->clear();
    hitt->clear();
    hitc->clear();
    tubeIds->clear();
    hitsType->clear();
    vertex_x->clear();
    vertex_y->clear();
    vertex_z->clear();
    vertex_t->clear();
    direction_x->clear();
    direction_y->clear();
    direction_z->clear();
    energy->clear();
  }

  outtree->Write();
  outfile->Close();

  return neventsSignal;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Get the hits information from a WCSim file object
 * @brief This code implement the "Sliding Window Algorithm" (<-> SWA) which extracts a trigger windows from the 1ms windows.
 *
 * The hits have to be time ordered, thus the use of the "sort" vector.
 * At this point we just extract windows, we do not care about the number of hits from signal or from noise.
 *
 * "Sliding Window Algorithm": roughtly, we "slide" a window of a given time (in ns) on the 1ms file, and if the number of hits contained in the window increases we keep sliding else we extract the window.
 * Thus, we extract non-overlapping windows that has the maximum hits.
 *
 * @param infilename Name of the file containing the events merged with the noise
 * @param outfilename Name of the file extracted using the "Sliding Window Algorithm"
 */

void get_hits_from_merged_file( string infilename, string outfilename, 
                                vector<TBox*> *vec_box, vector<double> *vec_BG, vector<double> *vec_SIG, vector<double> *vec_VTX,
                                int npmts,  float tankradius, float tankhalfz, float darkrate, int minimum_hit_signal, int nhit_per_200ns, int max_event_sig, int max_event_bg)
{
  int nPatch=0, stop_at_patch=0;  // Number of patch (1 patch = large noise event with short signal event pasted over it) and where we end the window extraction
  int NbSigInWindow = 0;          // Number of hits of signal in a trigger window
  int nWin_BG=0, nWin_SIG=0;      //< Number of events extracted with the SWA 
  int nWin=0, nWin_1st=0;         //< Number of extracted windows
  int shift;
  int windowType = -1;
  float t0=0., t1=0.;
  double dr, dx, dy, dz, dt, dR, dTh, dZ, dT, dE;
  double dirx, diry, dirz;
  bool IsNoiseOnly = true;
  bool ToExtract = false;
  vector<TimeWindow> tabWin;    //Table of index of the first hit ([0]), number of hits ([1]) in the windows, and if the type of the window ([3]: -1 for noise only, event # otherwise). Each slot can potentially hold a trigger window
  TRandom *r1 = new TRandom1();

  // Declare variables for branches of the input TTree
  int n_hits_IN;                        //< Number of hits in event
  vector<int> *eventNumber_IN = 0;     //< Event number used to optimize the SWA
  vector<float> *hitt_IN = 0;          //< Hits timings
  vector<float> *hitc_IN = 0;        //< Hits charges
  vector<float> *hitx_IN = 0;           //< x coordinates
  vector<float> *hity_IN = 0;           //< y coordinates
  vector<float> *hitz_IN = 0;           //< z coordinates
  vector<int> *tubeIds_IN = 0;
  vector<int> *hitsType_IN = 0;
  vector<float> *vertex_x_IN = 0;
  vector<float> *vertex_y_IN = 0;
  vector<float> *vertex_z_IN = 0;
  vector<float> *vertex_t_IN = 0;
  vector<float> *direction_x_IN = 0;
  vector<float> *direction_y_IN = 0;
  vector<float> *direction_z_IN = 0;
  vector<float> *energy_IN = 0;

  // Declare variables for branches of the output TTree
  int maxFrequency = 0;
  int mostFrequentEvent = 0;                                          //< Number of the event in the simulation
  int mostFrequentEvent_index = 0;                                    //< Index of the event with respect to the current patch
  int mostFrequentEvent_refIndex = 0;                                 //< Index in the simulation of the last event for each patch
  int mostFrequentEvent_refIndex_for_next_patch = 0;                  //< Index in the simulation of the first event for each patch
  int eventType;                                                      //< Type associated to the trigger window (noise=0, event=1)
  int n_hits;                                                         //< Number of hits in event
  float energy;                                                       //< Number of hits in event
  float vertex_x, vertex_y, vertex_z, vertex_t;                       //< x, y, z, t coordinates of gamma vertex
  float direction_x, direction_y, direction_z;                        //< Particle track start and stop position, and direction (start to stop)
  float hitx[MAX_HITS], hity[MAX_HITS], hitz[MAX_HITS], hitc[MAX_HITS], hitt[MAX_HITS];   //< x, y, z, charge and time arrays of hits coordinates
  int tubeIds[MAX_HITS];                                          //< Number of hits in event
  int hitsType[MAX_HITS];                                         //< Takes only 0 or 1 as a value, to indicate an actual value (1) or noise (0)
  int eventNumber[MAX_HITS];                                     //< Takes -1 to #signal as a value, to indicate an actual value (>=0) or noise (-1)


  TFile *file_IN = TFile::Open(infilename.data(), "READ");
  TTree *tree_IN = (TTree*)file_IN->Get("THits");
         tree_IN->SetBranchAddress("vertex_x", &vertex_x_IN);
         tree_IN->SetBranchAddress("vertex_y", &vertex_y_IN);
         tree_IN->SetBranchAddress("vertex_z", &vertex_z_IN);
         tree_IN->SetBranchAddress("vertex_t", &vertex_t_IN);
         tree_IN->SetBranchAddress("direction_x", &direction_x_IN);
         tree_IN->SetBranchAddress("direction_y", &direction_y_IN);
         tree_IN->SetBranchAddress("direction_z", &direction_z_IN);
         tree_IN->SetBranchAddress("energy", &energy_IN);
         tree_IN->SetBranchAddress("n_hits", &n_hits_IN);
         tree_IN->SetBranchAddress("eventNumber", &eventNumber_IN);
         tree_IN->SetBranchAddress("hitsType", &hitsType_IN);
         tree_IN->SetBranchAddress("hitx", &hitx_IN);
         tree_IN->SetBranchAddress("hity", &hity_IN);
         tree_IN->SetBranchAddress("hitz", &hitz_IN);
         tree_IN->SetBranchAddress("hitt", &hitt_IN);
         tree_IN->SetBranchAddress("hitc", &hitc_IN);
         tree_IN->SetBranchAddress("tubeIds", &tubeIds_IN);

  // Prepare output file and output tree - 1 outfile only with the same data
  TFile *outfile = new TFile(outfilename.data(), "recreate");
  TTree *outtree = new TTree("THits", "PMT hits from event + background splitted in windows");
         outtree->Branch("eventType", &eventType, "eventType/I");
         outtree->Branch("vertex_x", &vertex_x, "vertex_x/F");
         outtree->Branch("vertex_y", &vertex_y, "vertex_y/F");
         outtree->Branch("vertex_z", &vertex_z, "vertex_z/F");
         outtree->Branch("vertex_t", &vertex_t, "vertex_t/F");
         outtree->Branch("direction_x", &direction_x, "direction_x/F");
         outtree->Branch("direction_y", &direction_y, "direction_y/F");
         outtree->Branch("direction_z", &direction_z, "direction_z/F");
         outtree->Branch("energy", &energy, "energy/F");
         outtree->Branch("n_hits", &n_hits, "n_hits/I");
         outtree->Branch("eventNumber", eventNumber, "eventNumber[n_hits]/I");
         outtree->Branch("hitsType", hitsType, "hitsType[n_hits]/I");
         outtree->Branch("hitx", hitx, "hitx[n_hits]/F");
         outtree->Branch("hity", hity, "hity[n_hits]/F");
         outtree->Branch("hitz", hitz, "hitz[n_hits]/F");
         outtree->Branch("hitt", hitt, "hitt[n_hits]/F");
         outtree->Branch("hitc", hitc, "hitc[n_hits]/F");
         outtree->Branch("tubeIds", tubeIds, "tubeIds[n_hits]/I");


/*************************************************************************/
//Store hit info of the first patch for the plot

  tree_IN->GetEntry(0);

  //Fill time of hits of noise and of signal
  for (int j=0 ; j<n_hits_IN ; j++)
  {
    if(hitt_IN->at(j) > DRAW_TIME){break;}                                 //Stop when reaching time limit for drawing
    else
    {
      if((*eventNumber_IN)[j] == -1){vec_BG->push_back((*hitt_IN)[j]);}     // Fill the histogram from the vector
      else                          {vec_SIG->push_back((*hitt_IN)[j]);}    // Fill the histogram from the vector
    }
  }

  //Fill vertex time
  for (int j=0 ; j<vertex_t_IN->size() ; j++)
  {
    if(vertex_t_IN->at(j) > DRAW_TIME){break;}                                 //Stop when reaching time limit for drawing
    else{vec_VTX->push_back((*vertex_t_IN)[j]);}
  }


/*************************************************************************/
// Main loop over all events using the Sliding Window Algorithm

  cout << "\n\n" << endl;

  nPatch = tree_IN->GetEntries();
  stop_at_patch = nPatch;
  for (int i=0 ; i<nPatch ; i++)
  {
#ifdef VERBOSE
    cout << "\nProcessing patch #" << i << endl;
#endif

    tree_IN->GetEntry(i);
    if(nPatch >= 1){mostFrequentEvent_refIndex = mostFrequentEvent_refIndex_for_next_patch;}    //Set to the number of event that have passed until this patch
    mostFrequentEvent_refIndex_for_next_patch += energy_IN->size();                             //Increase by the number of signal events in the patch BUT this will used only at the next patch

    //Sliding Window Algorithm used to simulate triggers over the window of events
    vector<TimeWindow> tabWin_tmp = AlgTrigger(n_hits_IN, hitt_IN, nhit_per_200ns, npmts, darkrate);
    nWin = tabWin_tmp.size();
    if(i == 0){nWin_1st = nWin;}

    //Check if a window is a "signal" or a "noise" window
    for (int i_window=0 ; i_window<nWin ; i_window++)
    {
      shift  = tabWin_tmp[i_window].Index;          //Index of the first hit in the window
      n_hits = tabWin_tmp[i_window].Hits;           //Number of hits in the window
      // t0     = (*hitt_IN)[shift];                //Starting time of the window
      t0     = tabWin_tmp[i_window].TriggerTime;    //Trigger time for the window
      IsNoiseOnly = true;                           //False when a signal is identified 
      ToExtract = false;
      eventType = 0;
      NbSigInWindow = 0;

      //Check the number of hit of signal in the current window
      for (int j=0 ; j<n_hits ; j++)
      {
        hitt[j]        = (*hitt_IN)[shift + j] - t0;
        hitx[j]        = (*hitx_IN)[shift + j];
        hity[j]        = (*hity_IN)[shift + j];
        hitz[j]        = (*hitz_IN)[shift + j];
        hitc[j]        = (*hitc_IN)[shift + j];
        hitsType[j]    = (*hitsType_IN)[shift + j];
        eventNumber[j] = (*eventNumber_IN)[shift + j];

        if (hitsType[j]==1)
        {
          NbSigInWindow++;
          if(NbSigInWindow >= minimum_hit_signal){IsNoiseOnly = false;}     //Trigger condition
          // if(1.*NbSigInWindow >= 0.1*n_hits){IsNoiseOnly = false;}     //Mimic AB
        }
      }

      //No significant signal found after checking all the hits of a trigger window: then it's a "only noise" window
      if (IsNoiseOnly==true  &&  nWin_BG<max_event_bg)
      {
        ToExtract = true;
        nWin_BG++;
        eventType = 0;
        tabWin_tmp[i_window].Type = -1;

        //Fill the vertex branches with seemingly random values
        //NEED TO BE REPLACED BY RECONSTRUCTED VARIABLES !!!
        dR = r1->Uniform(0, 1) * tankradius;
        dTh= r1->Uniform(0, 1) * 2.*M_PI;
        dZ = r1->Uniform(-1, 1) * tankhalfz;
        dT = abs(r1->Gaus(0, 100));
        dE = abs(r1->Gaus(5, 3));
        r1->Sphere(dirx, diry, dirz, 1.);
        vertex_x   = dR * cos(dTh);
        vertex_y   = dR * sin(dTh);
        vertex_z   = dZ;
        vertex_t   = dT;
        direction_x= dirx;
        direction_y= diry;
        direction_z= dirz;
        energy     = dE;

        // cout << n_hits << endl;
      }

      //Significant signal found !!!
      else if (IsNoiseOnly==false  &&  nWin_SIG<max_event_sig)
      {
        ToExtract = true;
        nWin_SIG++;

        //Now we associate to the window the event generating the most hits in the window
        eventType = 1;
        maxFrequency = 0;
        mostFrequentEvent = -1;
        unordered_map<int, int> frequencyMap;                                  //Create an unordered_map to store the frequency of each signal number
        for (int j=0 ; j<n_hits ; j++)
        {
          if(eventNumber[j] >= 0)  //Noise hits have eventNumber==-1
          {
            frequencyMap[eventNumber[j]]++;       //Count the frequency of each event based on number of hit (that are not noise hits) in the window
          }
        }
        for (const auto& pair : frequencyMap) 
        {
            if (pair.first != -1   &&   pair.second > maxFrequency)
            {
                maxFrequency      = pair.second;
                mostFrequentEvent = pair.first;
            }
        }


        //Fill the vertex branches based on the info from the most frequent event and ONYL IF THERE IS A SIGNAL EVENT 
        if(mostFrequentEvent == -1)
        {
          ToExtract = false;    // Here, this is actually a noise event so we don't extract it
        }
        else
        {
          tabWin_tmp[i_window].Type = mostFrequentEvent;

          mostFrequentEvent_index = mostFrequentEvent - mostFrequentEvent_refIndex;

          vertex_x   = vertex_x_IN->at(mostFrequentEvent_index);
          vertex_y   = vertex_y_IN->at(mostFrequentEvent_index);
          vertex_z   = vertex_z_IN->at(mostFrequentEvent_index);
          vertex_t   = vertex_t_IN->at(mostFrequentEvent_index) - t0;
          direction_x= direction_x_IN->at(mostFrequentEvent_index);
          direction_y= direction_y_IN->at(mostFrequentEvent_index);
          direction_z= direction_z_IN->at(mostFrequentEvent_index);
          energy     = energy_IN->at(mostFrequentEvent_index);

          // cout << n_hits << "  " << NbSigInWindow << endl;
        }
      }


      //Fill data in TTree
      if (ToExtract == true)
      {
        outtree->Fill();
        tabWin.push_back(tabWin_tmp[i_window]);

#ifdef VERBOSE
        if ((nWin_BG + nWin_SIG) % PRINT_EVERY == 0)
        {
          cout << "Extracted trigger window #" << nWin_BG + nWin_SIG     //+1 to start counting at 1
               << ": [" << (*hitt_IN)[tabWin_tmp[i_window].Index] 
               << ", "     << (*hitt_IN)[tabWin_tmp[i_window].Index + tabWin_tmp[i_window].Hits - 1] 
               << "] ns, trigger time at " << tabWin_tmp[i_window].TriggerTime
               << " ns, " << tabWin_tmp[i_window].Hits
               << " hits,  event #"  << tabWin_tmp[i_window].Type << " (index " << mostFrequentEvent_index << " in current patch)";
               if(tabWin_tmp[i_window].Type >= 0){cout << "  with  Total energy="  << energy << " MeV";}
          cout << endl;
        }
#endif
      }
    }


    //Max number of graph reached, go out of loop
    if(nWin_BG>=max_event_bg  &&  nWin_SIG>=max_event_sig){stop_at_patch=i; break;}
  }


  cout << "Sliding Window Algorithm ending at patch #" << stop_at_patch << " after triggering " << nWin_SIG << " signal windows and " << nWin_BG << " noise windows" << endl;



  // Save and close out file
  outtree->Write();
  outfile->Close();



/*************************************************************************/
//Prepare histo output

  if(nWin_1st > nWin_SIG + nWin_BG){nWin_1st = nWin_SIG + nWin_BG;}
  cout << "\nPloting the first " << nWin_1st << " trigger windows over the timeline" << endl;

  //Plots only the windows for the first patch
  tree_IN->GetEntry(0);
  for (int i_window=0 ; i_window<nWin_1st ; i_window++)
  {
    shift  = tabWin[i_window].Index;          //Index of the first hit in the window
    n_hits = tabWin[i_window].Hits;           //Number of hits in the window
    windowType = tabWin[i_window].Type;       //Type of the window (noise or event)
    t0     = (*hitt_IN)[shift];               //Time of the first hit of the window
    t1     = (*hitt_IN)[shift + n_hits - 1];  //Time of the last hit of the window

    //Stop when reaching time limit for drawing
    if(t1 > DRAW_TIME){break;}

    // Create a light-colored rectangle
    else
    {
      TBox* rect = new TBox(t0, 0., t1, DRAW_HITY);    //Assumes vertical range from 0 to 1
      if(windowType == -1){rect->SetFillColor(kYellow);}            //Noise type
      else                {rect->SetFillColor(kCyan-9);}
      vec_box->push_back(rect);
    }
  }
}















// void get_hits_from_merged_file_vAB(string infilename, string outfilename, 
//                                 vector<TBox*> *vec_box, vector<double> *vec_BG, vector<double> *vec_SIG, vector<double> *vec_VTX, 
//                                 int npmts, float tankradius, float tankhalfz, float darkrate, int minimum_hit_signal, int nhit_per_200ns, int max_event)
// {
//   int nPatch = 0;               // Number of patch (1 patch = large noise event with short signal event pasted over it)
//   int NbSigInWindow = 0;        // Number of hits of signal in a trigger window
//   int nWin_BG=0, nWin_SIG=0;    //< Number of events extracted with the SWA 
//   int nWin=0, nWin_1st=100;     //< Number of extracted windows
//   int shift;
//   int windowType = -1;
//   float t0=0., t1=0.;
//   double dr, dx, dy, dz, dR, dTh, dZ;
//   bool IsNoiseOnly = true;
//   bool ToExtract = false;
//   vector<TimeWindow> tabWin;    //Table of index of the first hit ([0]), number of hits ([1]) in the windows, and if the type of the window ([3]: -1 for noise only, event # otherwise). Each slot can potentially hold a trigger window
//   TRandom *r1  = new TRandom1();
//   TRandom3 *r3 = new TRandom3();

//   // Declare variables for branches of the input TTree
//   Int_t n_hits_IN;                        //< Number of hits in event
//   float hitt_IN[MAX_HITS];          //< Hits timings
//   float hitc_IN[MAX_HITS];        //< Hits charges
//   float hitx_IN[MAX_HITS];           //< x coordinates
//   float hity_IN[MAX_HITS];           //< y coordinates
//   float hitz_IN[MAX_HITS];           //< z coordinates
//   int hitsType_IN[MAX_HITS];
//   float vertex_IN[3];
//   float vertex_0_IN[3];

//   // Declare variables for branches of the output TTree
//   int eventType;                                                      //< Type associated to the trigger window (noise=0, event=1)
//   int n_hits;                                                         //< Number of hits in event
//   float vertex[3];                                                    //< x, y, z coordinates of gamma vertex
//   float vertex_0[3];                                                  //< x, y, z coordinates of randomly reconstructed neutron vertex
//   float hitx[MAX_HITS], hity[MAX_HITS], hitz[MAX_HITS];   //< x, y, z arrays of hits coordinates
//   float hitc[MAX_HITS], hitt[MAX_HITS];                     //< Charge and time values of each hit
//   int hitsType[MAX_HITS];                                         //< Takes only 0 or 1 as a value, to indicate an actual value (1) or noise (0)


//   TFile *file_IN = TFile::Open(infilename.data(), "READ");
//   TTree *tree_IN = (TTree*)file_IN->Get("sk2p2");
//          tree_IN->SetBranchAddress("hitx", &hitx_IN);
//          tree_IN->SetBranchAddress("hity", &hity_IN);
//          tree_IN->SetBranchAddress("hitz", &hitz_IN);
//          tree_IN->SetBranchAddress("hitc", &hitc_IN);
//          tree_IN->SetBranchAddress("hitt",   &hitt_IN);
//          tree_IN->SetBranchAddress("n_hits", &n_hits_IN);
//          tree_IN->SetBranchAddress("hitsType", &hitsType_IN);
//          tree_IN->SetBranchAddress("vertex", &vertex_IN);
//          tree_IN->SetBranchAddress("vertex_0", &vertex_0_IN);

//   // Prepare output file and output tree - 1 outfile only with the same data
//   TFile *outfile = new TFile(outfilename.data(), "recreate");
//   TTree *outtree = new TTree("signal_and_background_PMT_hits_windows", "PMT hits from event + background splitted in windows");
//          outtree->Branch("eventType", &eventType, "eventType/I");
//          outtree->Branch("hitx", hitx, "hitx[n_hits]/F");
//          outtree->Branch("hity", hity, "hity[n_hits]/F");
//          outtree->Branch("hitz", hitz, "hitz[n_hits]/F");
//          outtree->Branch("hitc", hitc, "hitc[n_hits]/F");
//          outtree->Branch("hitt", hitt, "hitt[n_hits]/F");
//          outtree->Branch("n_hits", &n_hits, "n_hits/I");
//          outtree->Branch("hitsType", hitsType, "hitsType[n_hits]/I");
//          outtree->Branch("vertex", vertex, "vertex[3]/F");
//          outtree->Branch("vertex_0", vertex_0, "vertex_0[3]/F");


// /*************************************************************************/
// //Store hit info of the first patch for the plot


//   //Fill time of hits of noise and of signal
//   for (int i=0 ; i<nWin_1st ; i++)
//   {
//     tree_IN->GetEntry(i);
//     bool IsSignal = false;

//     for (int j=0 ; j<n_hits_IN ; j++)
//     {
//       if(hitsType_IN[j] == 0){vec_BG->push_back(hitt_IN[j] + i*500 + 10);}     // Fill the histogram from the vector
//       else                   {vec_SIG->push_back(hitt_IN[j] + i*500 + 10);  IsSignal = true;}    // Fill the histogram from the vector
//     }

//     if(IsSignal == true){vec_VTX->push_back(i*500 + 10);}
//   }



// /*************************************************************************/
// // Main loop over all events using the Sliding Window Algorithm

//   cout << "\n\n" << endl;

//   nWin = tree_IN->GetEntries();
//   int index_signal = 0;
//   int trigger_hit = 11;
//   int index_vtx = 0;

//   vector<float> check_if_event_is_unique = {0.};

//   //Check if a window is a "signal" or a "noise" window
//   // for (int i_window=0 ; i_window<nWin ; i_window++)
//   for (int i_window=10000000 ; i_window<nWin ; i_window++)
//   {
//     tree_IN->GetEntry(i_window);

//     TimeWindow window_tmp;
//               window_tmp.Index = 0;
//               window_tmp.Hits  = 0;
//               window_tmp.Type  = 0;
//               window_tmp.TriggerTime = 0.;

//     shift  = 0;      //Index of the first hit in the window
//     n_hits = n_hits_IN;       //Number of hits in the window
//     t0     = hitt_IN[shift + trigger_hit];          //Trigger time for the window
//     IsNoiseOnly = true;                   //False when a signal is identified (this is where the "trigger" selection should apply)
//     ToExtract = false;
//     eventType = 0;
//     NbSigInWindow = 0;

//     if( n_hits != 0   &&   n_hits < 70)
//     {
//       window_tmp.Index = shift;
//       window_tmp.Hits  = n_hits;
//       window_tmp.TriggerTime = t0;

//       //Check the number of hit of signal in the current window
//       for (int j=0 ; j<n_hits ; j++)
//       {
//         hitx[j]     = hitx_IN[shift + j];
//         hity[j]     = hity_IN[shift + j];
//         hitz[j]     = hitz_IN[shift + j];
//         hitc[j]     = hitc_IN[shift + j];
//         hitt[j]     = hitt_IN[shift + j];
//         hitsType[j] = hitsType_IN[shift + j];

//         // if (IsNoiseOnly  &&  hitsType[j]==1)
//         if (hitsType[j]==1)
//         {
//           NbSigInWindow++;
//           // if(NbSigInWindow >= minimum_hit_signal){IsNoiseOnly = false;}     //Trigger condition
//           if(1.*NbSigInWindow > 0.1*n_hits){IsNoiseOnly = false;}     //Mimic AB
//         }
//       }

//       //No significant signal found after checking all the hits of a trigger window: then it's a "only noise" window
//       if (IsNoiseOnly==true  &&  nWin_BG<max_event)
//       {
//         dR = sqrt(r3->Uniform()) * tankradius;
//         dTh= r3->Uniform() * 2.*M_PI;
//         dZ = (r3->Uniform()*2.-1.) * tankhalfz;
//         dR = r1->Gaus(0, 50);
//         r1->Sphere(dx, dy, dz, dr);
//         float value_to_check = dR * cos(dTh) + dx;
//         auto it = std::find(check_if_event_is_unique.begin(), check_if_event_is_unique.end(), value_to_check);
//         if (it == check_if_event_is_unique.end())
//         {
//           check_if_event_is_unique.push_back(value_to_check);    // Value doesn't exist, add it to the vector

//           ToExtract = true;
//           nWin_BG++;
//           eventType = 0;
//           window_tmp.Type = -1;

//           //Fill the vertex branches with super negative value to identify it
//           vertex[0]   = dR * cos(dTh);
//           vertex[1]   = dR * sin(dTh);
//           vertex[2]   = dZ;
//           vertex_0[0] = dR * cos(dTh) + dx;
//           vertex_0[1] = dR * sin(dTh) + dy;
//           vertex_0[2] = dZ + dz;
//           // cout << n_hits << endl;
//         }
//       }

//       //Significant signal found !!!
//       else if(IsNoiseOnly==false  &&  nWin_SIG<max_event)
//       {
//         // Check if the neutron vertex associated to the signal has already been used
//         float value_to_check = vertex_0_IN[0];
//         auto it = std::find(check_if_event_is_unique.begin(), check_if_event_is_unique.end(), value_to_check);
//         if (it == check_if_event_is_unique.end())
//         {
//           check_if_event_is_unique.push_back(value_to_check);    // Value doesn't exist, add it to the vector

//           ToExtract = true;
//           nWin_SIG++;
//           eventType = 1;
//           window_tmp.Type = index_vtx;

//           //Fill the vertex branches based on the indo from the most frequent event
//           vertex[0]   = vertex_IN[0];
//           vertex[1]   = vertex_IN[1];
//           vertex[2]   = vertex_IN[2];
//           vertex_0[0] = vertex_0_IN[0];
//           vertex_0[1] = vertex_0_IN[1];
//           vertex_0[2] = vertex_0_IN[2];
//           index_vtx++;

//           // cout << n_hits << "  " << NbSigInWindow << endl;
//         }
//       }

//       else
//       {
//         ToExtract   = false;
//         vertex[0]   = -99999.;
//         vertex[1]   = -99999.;
//         vertex[2]   = -99999.;
//         vertex_0[0] = -99999.;
//         vertex_0[1] = -99999.;
//         vertex_0[2] = -99999.;
//       }


//       //Fill data in TTree
//       if (ToExtract == true)
//       {
//         outtree->Fill();
//         tabWin.push_back(window_tmp);

// #ifdef VERBOSE
//         if ((nWin_BG + nWin_SIG) % PRINT_EVERY == 0)
//         {
//           cout << "Stored trigger window #" << nWin_BG + nWin_SIG     //+1 to start counting at 1
//                << " at event #" << i_window     //+1 to start counting at 1
//                << ": [" << hitt_IN[window_tmp.Index] 
//                << ", "     << hitt_IN[window_tmp.Index + window_tmp.Hits - 1] 
//                << "] ns, trigger time at " << window_tmp.TriggerTime
//                << " ns, " << window_tmp.Hits
//                << " hits,  event type="  << window_tmp.Type;
//           cout << endl;
//         }
// #endif
//       }

//       if(nWin_SIG>=max_event  &&  nWin_BG>=max_event)
//       {
//         cout << "Ending after having reprocessed " << nWin_SIG << " signal windows and " << nWin_BG << " noise windows" << endl;
//         break;
//       }
//     }
//   }


//   // Save and close out file
//   outtree->Write();
//   outfile->Close();


// /*************************************************************************/
// //Prepare histo output

//   if(nWin_1st > nWin_SIG + nWin_BG){nWin_1st = nWin_SIG + nWin_BG;}
//   cout << "\nPloting the first " << nWin_1st << " trigger windows over the timeline" << endl;

//   //Plots only the windows for the first patch
//   for (int i_window=0 ; i_window<nWin_1st ; i_window++)
//   {
//     tree_IN->GetEntry(i_window);

//     shift  = tabWin[i_window].Index;          //Index of the first hit in the window
//     n_hits = tabWin[i_window].Hits;           //Number of hits in the window
//     windowType = tabWin[i_window].Type;       //Type of the window (noise or event)
//     t0     = hitt_IN[shift] + i_window*500 + 10;                 //Time of the first hit of the window
//     t1     = hitt_IN[shift + n_hits - 1] + i_window*500 + 10;    //Time of the last hit of the window

//     //Stop when reaching time limit for drawing
//     if(t1 > DRAW_TIME){break;}

//     // Create a light-colored rectangle
//     else
//     {
//       TBox* rect = new TBox(t0, 0., t1, DRAW_HITY);    //Assumes vertical range from 0 to 1
//       if(windowType == -1){rect->SetFillColor(kYellow);}            //Noise type
//       else                {rect->SetFillColor(kCyan-9);}
//       // rect->SetFillStyle(3004);
//       vec_box->push_back(rect);
//     }
//   }
// }




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv) 
{
  int index_box_signal=-1, index_box_noise=-1;
  int hSize = 500;                //Number of bins to split the time intervalle into
  int max_event_sig = -1;
  double time_gap = TIME_GAP;     //Average time gap between signal, in ns
  vector<TBox*> vec_box;
  vector<double> vec_BG;
  vector<double> vec_SIG;
  vector<double> vec_VTX;

  string path="", radical_signal="", radical_noise="";
  string input_signal="", output_signal="", input_noise="", output_noise="", output_merged="", output_windows="";
  string canvas_name = "/sps/t2k/lperisse/Soft/Preprocess_WCSim/plots/SWA_event_selection_scheme.png";

  //Fetch arguments
  cout << "\n" << endl;
  arguments arglist = FetchInput(argc, argv);
  string name_noise = arglist.inputNoise;
  string name_signal = arglist.inputSignal;
  if(arglist.signalRate){time_gap = TIME_GAP;}
  int event_type = arglist.eventType; 
  int nb_PMT_ID = arglist.PMT;
  float tankradius = arglist.tankRadius;
  float tankhalfz  = arglist.tankHalfz;
  float darkrate   = arglist.darkRate;
  int minimum_hit_signal = arglist.windowSignal;
  int nhit_per_200ns = arglist.nhits200ns;
  int max_event_bg = arglist.maxEventBg;
  // arglist.verbose;
  cout << "\n" << endl;

  //Set all input and output file names
  if(name_signal != "")
  {
    radical_signal = name_signal.substr(0, name_signal.size()-5);         //name of the file without the ".root"
    input_signal   = Form("%s",           name_signal.c_str());           //Signal files
    output_signal  = Form("%s_POST.root", radical_signal.c_str());        //Signal files once analyzed
  }
  if(name_noise  != "")
  {
    radical_noise  = name_noise.substr(0, name_noise.size()-5);           //name of the file without the ".root"
    input_noise    = Form("%s",           name_noise.c_str());            //Noise files
    output_noise   = Form("%s_POST.root", radical_noise.c_str());         //Noise files once analyzed
  }
  if(name_signal != ""   &&   name_noise != "")
  {
    output_merged  = Form("%s_MERGED.root",  radical_signal.c_str());     //Merged signal + noise file
    output_windows = Form("%s_WINDOWS.root", radical_signal.c_str());     //Windowed signal + noise events of noise only events
  }

/*************************************************************************/
//Post-simulation steps

  // if(name_signal != ""){get_hits_from_WCSim(input_signal, output_signal, event_type);}     //Select the events hits and vertex info
  // if(name_noise  != ""){get_hits_from_WCSim(input_noise, output_noise, 0);}       //Select noise hits
  if(name_signal != ""   &&   name_noise != ""){max_event_sig = merge_signal_to_noise(output_signal, output_noise, output_merged, time_gap);}                //Merge event and noise dataset by inserting events on top of noise window
  if(max_event_sig > 0){get_hits_from_merged_file(output_merged, output_windows, &vec_box, &vec_BG, &vec_SIG, &vec_VTX, nb_PMT_ID, tankradius, tankhalfz, darkrate, minimum_hit_signal, nhit_per_200ns, max_event_sig, max_event_bg);}     //Extact windows containing either event + noise events or noise only events, using a sliding window algorithm
  else{cout << "BEWARE, NO SIGNAL EVENT CAN BE EXTRACTED FROM THE MERGED FILE !!!" << endl;}



/*************************************************************************/
//Specific lines to reprocess Antoine Beauchene files

  // output_merged  = "/sps/t2k/lperisse/skynet_Lorenzo/Datasets/dataset_SKDetSim_2M_photons/dataset.root";
  // output_windows = "/sps/t2k/lperisse/skynet_Lorenzo/Datasets/WCSIM_outputs/dataset_SKDetSim_2M_photons_WINDOWS_bench_LP_fix2.root";
  // output_windows = "/sps/t2k/lperisse/skynet_Lorenzo/Datasets/WCSIM_outputs/dataset_SKDetSim_2M_photons_WINDOWS_bench_LP_eval.root";
  // output_merged  = "/sps/t2k/lperisse/skynet_Lorenzo/Datasets/dataset_train_SKDetSim_more_info/dataset.root";
  // output_windows = "/sps/t2k/lperisse/skynet_Lorenzo/Datasets/WCSIM_outputs/dataset_train_SKDetSim_more_info_WINDOWS_bench_LP_fix.root";

  // get_hits_from_merged_file_vAB(output_merged, output_windows, &vec_box, &vec_BG, &vec_SIG, &vec_VTX, nb_PMT_ID, tankradius, tankhalfz, darkrate, minimum_hit_signal, nhit_per_200ns, max_signal);     //Extact windows containing either event + noise events or noise only events, using a sliding window algorithm



/*************************************************************************/
//Plot the SWA event selection scheme (=timeline of hits and events)

  //Get trigger threshold
  double time_min = vec_BG[0];
  double dark_rate_Hz = darkrate * 1000.;                                          //Dark rate of PMTs in Hz
  double trigger_window_seconds = TIME_WINDOW * 1E-9;                              //Length of the trigger window in seconds
  double average_occupancy = dark_rate_Hz * trigger_window_seconds * nb_PMT_ID;    //Average number of dark rate hits in a trigger window
  double size_window  = TIME_WINDOW;                                               //Size of the trigger window in ns
  double n_windows  =  (DRAW_TIME - time_min) / TIME_WINDOW;                       //Number of time windows that could fit into the total plotted time
  double nbin_per_window = (hSize / n_windows);                                    //Number of bin per window
  double ndigitsThreshold = 1.*round(average_occupancy  +  nhit_per_200ns*TIME_WINDOW/200.) / nbin_per_window;    //Number of hits in a trigger window to issue trigger (for instance for SK: ndigitsThreshold = 25 for size_window = 200)

  TH1F *hist_Blank= new TH1F("Blank",    "Scheme of SWA event selection", hSize, time_min, DRAW_TIME);      //Background noise only
  TH1F *hist_BG   = new TH1F("BG_hits",  "BG_hits",  hSize, time_min, DRAW_TIME);      //Background noise only
  TH1F *hist_SIG  = new TH1F("SIG_hits", "SIG_hits", hSize, time_min, DRAW_TIME);      //Signal superimposed to noise

  //Set background and signal histograms
  for(int i=0 ; i<vec_BG.size() ; i++){hist_BG->Fill(vec_BG[i]);    hist_SIG->Fill(vec_BG[i]);}  //Fill noise hit time
  for(int i=0 ; i<vec_SIG.size() ; i++){hist_SIG->Fill(vec_SIG[i]);}                             //Fill signal hit time

  double hLowEdge  = hist_Blank->GetBinLowEdge(1);
  double hHighEdge = hist_Blank->GetBinLowEdge(hist_Blank->GetNbinsX()+1);
  double mean_nHitsNoise = hist_BG->Integral() / hist_BG->GetNbinsX();  // Compute the mean of the number of hits in noise events
  double std_nHitsNoise  = .0;
    for(int i=0 ; i<hist_BG->GetNbinsX() ; i++){std_nHitsNoise += (hist_BG->GetBinContent(i+1) - mean_nHitsNoise)*(hist_BG->GetBinContent(i+1) - mean_nHitsNoise);}
  std_nHitsNoise = sqrt(std_nHitsNoise / hist_BG->GetNbinsX());


  TCanvas *Can = new TCanvas("Can_timeline_hits", "Can_timeline_hits", 1500, 500);
          Can->SetFillColor(10);
          // Can->SetGrid();
          gStyle->SetOptStat(0);

  //Draw signal histogram
  hist_Blank->SetLineColor(kWhite);
  hist_Blank->SetLineWidth(0);
  hist_Blank->GetXaxis()->SetRangeUser(hLowEdge, hHighEdge);
  hist_Blank->GetXaxis()->SetTitle("Time of hits  [ns]");
  hist_Blank->GetXaxis()->CenterTitle();
  hist_Blank->GetYaxis()->SetRangeUser(0., DRAW_HITY);
  hist_Blank->GetYaxis()->SetTitle("Number of hits");
  hist_Blank->GetYaxis()->CenterTitle();
  hist_Blank->Draw("hist");


  // //Draw all windowed events + last box = standard deviation of number of hits of noise
  for (int i=0 ; i<vec_box.size() ; i++)
  {
    vec_box[i]->Draw("same");
    if(index_box_signal==-1  &&  vec_box[i]->GetFillColor()==kCyan-9){index_box_signal = i;}
    if(index_box_noise==-1   &&  vec_box[i]->GetFillColor()==kYellow){index_box_noise  = i;}
  }

  //Draw vertex time
  TLine *line_VTX[vec_VTX.size()];
  for(int i=0 ; i<vec_VTX.size() ; i++)
  {
    line_VTX[i] = new TLine(vec_VTX[i], 0., vec_VTX[i], DRAW_HITY);
    line_VTX[i]->SetLineWidth(2);
    line_VTX[i]->SetLineColor(kBlue);
    line_VTX[i]->Draw("same");
  }

  //Draw superimposed signal noise histogram
  hist_SIG->SetLineColor(kRed);
  hist_SIG->SetLineWidth(2);
  hist_SIG->Draw("hist same");

  //Draw noise histogram
  hist_BG->SetLineColor(kGray);
  hist_BG->SetLineWidth(2);
  hist_BG->Draw("hist same");

  // Create a colored line for mean and standard deviations
  TLine* line_mean = new TLine(hLowEdge, mean_nHitsNoise, hHighEdge, mean_nHitsNoise);
         line_mean->SetLineWidth(1);
         line_mean->SetLineColor(kBlack);
         line_mean->Draw("same");

  TLine* line_stdDown = new TLine(hLowEdge, mean_nHitsNoise-std_nHitsNoise, hHighEdge, mean_nHitsNoise-std_nHitsNoise);
         line_stdDown->SetLineWidth(1);
         line_stdDown->SetLineColor(kBlack);
         line_stdDown->SetLineStyle(7);
         line_stdDown->Draw("same");

  TLine* line_stdUp = new TLine(hLowEdge, mean_nHitsNoise+std_nHitsNoise, hHighEdge, mean_nHitsNoise+std_nHitsNoise);
         line_stdUp->SetLineWidth(1);
         line_stdUp->SetLineColor(kBlack);
         line_stdUp->SetLineStyle(7);
         line_stdUp->Draw("same");

  TLine* line_stdUp2 = new TLine(hLowEdge, mean_nHitsNoise+2*std_nHitsNoise, hHighEdge, mean_nHitsNoise+2*std_nHitsNoise);
         line_stdUp2->SetLineWidth(1);
         line_stdUp2->SetLineColor(kBlack);
         line_stdUp2->SetLineStyle(2);
         line_stdUp2->Draw("same");

  TLine* line_trigger = new TLine(hLowEdge, ndigitsThreshold, hHighEdge, ndigitsThreshold);
         line_trigger->SetLineWidth(1);
         line_trigger->SetLineColor(kOrange);
        //  line_trigger->SetLineStyle(2);
         line_trigger->Draw("same");

  TLegend *Legend = new TLegend(0.75, 0.5, 0.9, 0.9, NULL,"brNDC");
           if(index_box_noise > -1){Legend->AddEntry(vec_box[index_box_noise],  "Window noise ",  "f");}
           Legend->AddEntry(hist_BG,     "Dark noise hits ");
           Legend->AddEntry(line_mean,   "Dark noise mean ",    "l");
           Legend->AddEntry(line_stdUp,  "Dark noise 1#sigma ", "l");
           Legend->AddEntry(line_stdUp2, "Dark noise 2#sigma ", "l");
           Legend->AddEntry(line_trigger,"Trigger threshold ", "l");
           Legend->AddEntry(vec_box[index_box_signal],  "Window signal ", "f");
           Legend->AddEntry(hist_SIG,    "Signal hits ");
           Legend->AddEntry(line_VTX[0], "Vertex time origin ", "l");
           Legend->SetTextSize(0.035);
           Legend->Draw();

  Can->Draw();
  Can->SaveAs(canvas_name.data());


  return 1;
}

