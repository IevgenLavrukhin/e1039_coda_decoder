#include <iostream>
#include <iomanip>

#include <TList.h>
#include <TH1.h>
#include <TF1.h>
#include <TString.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <stdio.h>
#include <stdlib.h>
#include <TObjString.h>

#include <map>
#include <vector>

#include "THaCodaFile.h"
#include "THaEtClient.h"

#define MAX_EVENT_SIZE 70000

using namespace std;

//Event storage
class Event
{
public:
    Event();
    void setEventID(int codaEventID, int eventID);
    void addHit(unsigned int hit);
    void addV1495Hit(std::vector <unsigned int> tdc_word, unsigned int common_stop);

    void setHeader(unsigned int header);
    void setV1495Header(unsigned int stop_time, unsigned int n_events);

    void fillInfo(unsigned int word, double triggerTime, int& channelID, double& tdcTime);
    double decodeTime(unsigned int word);

public:
    int codaEventID;
    int eventID;

    int nEntriesExp;
    double triggerTime;

    vector<unsigned int> channels;
    vector<double> tdcTimes;
};

//TDC storage
class TDC
{
public:
    TDC();
    void init();
    TString check();

    void finalizeEvent(int codaEventID, int eventID);
    void fillHeader(unsigned int header);
    void fillV1495Header(unsigned int stop_time, unsigned int n_events);
    void fillHit(unsigned int hit);
    void fillV1495Hit(std::vector <unsigned int> tdc_word, unsigned int common_stop);

public:
    int boardID;
    vector<Event> events;
};

//ROC storage
class ROC
{
public:
    void init();
    TString check();

public:
    int rocID;
    int nTDCs;
    vector<TDC> tdcs;
};

//===========================================================================================
const int NROCs = 15;

const int nV1495_Boards = 5;

//const int NROCs = 1;
unsigned int RocIDs[NROCs] = {12, 13, 14, 15, 17, 18, 19, 21, 22, 23, 25, 26, 28, 30, 31};
//unsigned int RocIDs[NROCs] = {6};
unsigned int NTDCs[NROCs]  = {6,  3,  5,  6,  7,  7,  6,  6,  6,  7,  nV1495_Boards,  7,  5,  7,  5 };
//unsigned int NTDCs[NROCs]  = {1};

// TDC mapper for v1495 TDCs:
unsigned int v1495_Board_ID[nV1495_Boards] = {0x0420, 0x0430, 0x0460, 0x0470, 0x480}; //L0_T, L0_B, L1_T, L1_B, L2


int get_v1495_number(unsigned int firmware_ID){
  for(int i=0; i<nV1495_Boards; i++){
    if(firmware_ID == v1495_Board_ID[i])
        return i;
  }

  return -1;
}




int main(int argc, char* argv[])
{
    THaCodaData* coda = new THaCodaFile(TString(argv[1]));

    //Initialization
    map<int, ROC> rocs;
    map<int, bool> ARMdead;
    bool ARMdeadFlag = false;
 	  vector<int> eventTys;
 	  for(int i = 0; i < NROCs; ++i)
    {
        ROC newROC;

        newROC.rocID = RocIDs[i];
        newROC.nTDCs = NTDCs[i];
        for(int j = 0; j < newROC.nTDCs; ++j)
        {
            TDC newTDC;
            newTDC.boardID = j;
            newROC.tdcs.push_back(newTDC);
        }

        rocs[RocIDs[i]] = newROC;
        rocs[RocIDs[i]].init();

        ARMdead[RocIDs[i]] = false;
    }

    //Book output tuple
    int rocID;
    int boardID;
    int eventID;
    int channelID;
    double tdcTime;
    int eventTy;

    TFile* saveFile = new TFile(argv[2], "recreate");
    TTree* saveTree = new TTree("save", "save");

    saveTree->Branch("rocID", &rocID);
    saveTree->Branch("boardID", &boardID);
    saveTree->Branch("channelID", &channelID);
    saveTree->Branch("eventID", &eventID);
    saveTree->Branch("tdcTime", &tdcTime);
    saveTree->Branch("eventTy", &eventTy);

    //Read & decode
    int spillID = 0;
    int targetPos = 0;
    int codaEventID = 1;
    int bosEventID = 0;
    int eosEventID = 0;
    int minSpillID = -1;


    int event_counter = 0;

    unsigned int* data;
    bool firstBOS = true;
    while(true)
    {
        event_counter ++;
        if(event_counter%100000 == 0){
            printf("Processing Event # : %i\n", event_counter);
        }

        int status = coda->codaRead();
        if(status != 0)
        {
            if(status == -1)
            {
                coda->codaClose();
                break;
            }
            else
            {
                cout << "Spotted a corruptted event." << endl;
                continue;
            }
        }

        data = coda->getEvBuffer();
        int eventType = data[1] >> 16;
        int nWordsTotal = data[0] + 1;

        if(eventType == 11 || eventType == 0x14) //BOS or normal end of run
        {
            //Run event check -- only when:
            //  1. spillID larger than minimum;
            //  2. not the first spill
            //  3. all ARM cores are working fine
            //cout << "spillID = "<< spillID<< ", eosEventID = "<<eosEventID<<", bosEventID = "<<bosEventID<<endl;
            if(spillID > minSpillID && !firstBOS && !ARMdeadFlag && eosEventID > bosEventID)
            {
				       // cout << "Spill " << spillID << "  BOS " << bosEventID << "  EOS " << eosEventID << "  targetPos " << targetPos << endl;
              //  for(int i = 0; i < NROCs; ++i)  cout << rocs[RocIDs[i]].check() << endl;  //print basic info

                //dump data to tuple
                unsigned int nEvents = rocs[RocIDs[0]].tdcs[0].events.size();
                for(unsigned int iEvt = 0; iEvt < nEvents; ++iEvt)
                {
                    for(unsigned int iRoc = 0; iRoc < NROCs; ++iRoc)
                    {
                        for(unsigned int iTDC = 0; iTDC < rocs[RocIDs[iRoc]].nTDCs; ++iTDC)
                        {
			                      if(iEvt >= rocs[RocIDs[iRoc]].tdcs[iTDC].events.size()) continue;

                            Event thisEvent = rocs[RocIDs[iRoc]].tdcs[iTDC].events[iEvt];

                            //printf("\n\n\n\n\n ROC = %i: nHits = %i, TDC = %i \n\n\n\n", RocIDs[iRoc], thisEvent.tdcTimes.size(), iTDC);

                            for(unsigned int iHit = 0; iHit < thisEvent.tdcTimes.size(); ++iHit)
                            {
                                rocID = RocIDs[iRoc];
                                boardID = rocs[rocID].tdcs[iTDC].boardID;
                                eventID = thisEvent.eventID;
                                channelID = thisEvent.channels[iHit];
                                tdcTime = thisEvent.tdcTimes[iHit];
								                eventTy = eventTys[iEvt]; //
                                saveTree->Fill();

                              //  printf("roc =%i, bord = %i, eventID = %i,  eventTY = %i \n", rocID, boardID, eventID, eventTy);
                            }
                        }
                    }
                }

                targetPos = 0;
				        eventTys.clear();
            }
            firstBOS = false;

            //Quit if ARM is dead
            if(ARMdeadFlag) return 1;
            //Clear storage and reset ARM status flag
            ARMdeadFlag = false;
            for(int i = 0; i < NROCs; ++i)
            {
                rocs[RocIDs[i]].init();
                ARMdead[RocIDs[i]] = false;
            }

            if(eventType == 11)
            {
                bosEventID = codaEventID;

                ++codaEventID;
                continue;
            }
            else
            {
                break;
            }
        }
        else if(eventType == 14){ //v1495 TDC data from FGPA trigger Roc15 => Roc25 here
//          printf("eventType == 14 => Size = %i \n", nWordsTotal);

          int tdc_id = -1;
          int ts_event_ID = -1;
          int triggerType = -1;
          int n_word = 7;

        //  cout << "Spill " << spillID << "  BOS " << bosEventID << "  EOS " << eosEventID  << endl;

          if(spillID > minSpillID)
          while (n_word < nWordsTotal){
            if(data[n_word] == 0x13378eef){
              tdc_id++;

              unsigned int b_ID = data[++n_word];
              unsigned int time_window = data[++n_word];
              unsigned int n_hits = data[++n_word] & 0xffff;
              unsigned int commot_stop = data[++n_word] & 0xfff;


              if(n_hits == 0xd1ad || commot_stop == 0xd2ad){ //if TDC srewed up readout it sends a garbage. Need to check it.
                n_hits  = 0;
                commot_stop = 0;
              }

              for (int i = 0; i<n_hits; i++){     //loop over TDC hits:
                unsigned int tdc_word = data[++n_word];
                int tdc_ch = (int) (tdc_word & 0xff00) >> 8;
                unsigned int tdc_time = commot_stop - (tdc_word & 0xff);

                if(tdc_ch > 95){
                //  printf("tdc id = %i => board ID = 0x%x; total # hits = %i => hit = %i ch = %i tdc = %i \n", tdc_id, b_ID, n_hits, i, tdc_ch, tdc_time);
                //  printf("ts_event_ID = %i, triggerTipe = %i \n",ts_event_ID, triggerType);
                }

              //fill the free:
                if(triggerType > 0 && ts_event_ID > 0){
                    rocID = 25; //RocIDs[iRoc];
              //      boardID = tdc_id;
              //      eventID = ts_event_ID;
              //      channelID = tdc_ch;
              //      tdcTime = (double) tdc_time;//*18.86/16.0; // need to conver it properly
              //      eventTy = triggerType; //
              //      saveTree->Fill();
                }

            } //end for loop
            n_word ++;

            }
            else if (data[n_word]==0xe906f00f){
               ts_event_ID = (int) data[++n_word];
               triggerType = (int)data[++n_word];
               n_word ++;
              // printf("ts_event_ID = %i, triggerTipe = %i \n",ts_event_ID, triggerType);

            }else{  //all other words are skipped ?
                n_word ++;

            }

          }  //end while loop

          ++codaEventID;

          continue;

        }
        else if(eventType == 129)   //spill counter
        {
            TString spillIDstr;
            for(int i = 4; i < nWordsTotal; ++i)
            {
                for(int j = 0; j < 4; ++j)
                {
                    spillIDstr = Form("%s%c", spillIDstr.Data(), (data[i] >> (j*8)) & 0xff);
                }
            }
            spillID = spillIDstr.Atoi();

            ++codaEventID;
            continue;
        }
        else if(eventType == 130) //Slow control
        {
            TString slowcontrolStr;
            for(int i = 4; i < nWordsTotal; ++i)
            {
                for(int j = 0; j < 4; ++j)
                {
                    slowcontrolStr = Form("%s%c", slowcontrolStr.Data(), (data[i] >> (j*8)) & 0xff);
                }
            }

            TObjArray* slowcontrolDataGroup = slowcontrolStr.Tokenize("\n");
            if(slowcontrolDataGroup->GetEntries() > 117)
            {
                TString targetString = ((TObjString*)(slowcontrolDataGroup->At(117)))->String();
                TObjArray* targetDataGroup = targetString.Tokenize(" ");
                if(targetDataGroup->GetEntries() == 4)
                {
                    targetPos = ((TObjString*)(targetDataGroup->At(2)))->String().Atoi();
                }
                delete targetDataGroup;
            }
            delete slowcontrolDataGroup;

            ++codaEventID;
            continue;
        }
        else if(eventType == 12 || eventType == 17 || eventType == 18 || eventType == 132 || eventType == 130 || eventType == 140)
        {
            if(eventType == 12) eosEventID = codaEventID;

            ++codaEventID;
            continue;
        }
        if(spillID <= minSpillID) continue;
         cout << " codaEventID = " << codaEventID << ", nWords = " << nWordsTotal << " " << eventType << endl;

         /*
         if(eventType ==10){
           for(int i=0; i<nWordsTotal; i++){
              if(i%10==0) printf("\n");
             printf("0x%x\t", data[i]);

           }
         }

         */

        int iWord = 7;
        while(iWord < nWordsTotal)
        {
            //entry per ROC
            int nWordsRoc = data[iWord++];
            int maxRocWordID = iWord + nWordsRoc;
            int rocID = (data[iWord++] & 0x00ff0000) >> 16;
            //cout << "RocID = " << dec << rocID << ", nWordsRoc = " << dec << nWordsRoc << endl;
            if(rocID == 25)// || rocID == 30)// || rocID == 2)
            {
                //cout << "RocID = " << dec << rocID << ", nWordsRoc = " << dec << nWordsRoc << endl;
            }

            ++iWord; ++iWord; ++iWord; //neglect the first 3 words
            while(iWord < maxRocWordID)
            {
                if(data[iWord] == 0xe906f00f) //trigger type from TS
                {
                    ++iWord;
					          unsigned int nEvents=0;
					          if (data[iWord]>0x00000000){
						             nEvents=(data[iWord]-1)/2;
						             //cout<<nEvents<<endl;
						        }
                    ++iWord;
					          for (unsigned int i=0; i<nEvents; i++){
						            eventTys.push_back(data[iWord]);
                    	  //cout<<data[iWord]<<" "<<data[iWord+1]<<endl;
						            ++iWord;
                      	++iWord;
					          }
                    //cout << data[iWord++] << "  " << data[iWord++] << endl;
                    //++iWord; ++iWord; //not needed for data check
                	iWord = maxRocWordID;
                }
                else if(data[iWord] == 0xe906f005) // V1495 TRigger TDC readout!!!!!
                {
                    unsigned int v1495_TDC_ID = data[++iWord]; //first word is TDC ID of the board
                    unsigned int n_v1495_TDC_words = data[++iWord]; //second word is number of word

                    int v1495_board_num = get_v1495_number(v1495_TDC_ID);

                    if(n_v1495_TDC_words != 0){
                    //  printf("EventTY: %i;\t V1495 TDC bank started for 0x%x => board number = %i =>  number of words = %i\n", \
                                          eventType, v1495_TDC_ID, v1495_board_num, n_v1495_TDC_words);


                      std::vector <unsigned int> v1495_hits;

                      int v1495extraWords=0; // this is needed to take into account 2 extra words per physics event (stop time & codaID)
                      int i=0;

                      std::vector<unsigned int> v1495_tdc;

                      while (i< n_v1495_TDC_words+v1495extraWords){ //up to 6 events per readout  & 2 extra words stop time & coda event ID
                        ++iWord;

                        //printf("data[%i] = 0x%x \n",i, data[iWord]);

                        if(data[iWord]>>16 == 0){
                        //  printf("TDC data[%i] = 0x%x \n",i, data[iWord]);
                          v1495_tdc.emplace_back(data[iWord]);
                        }

                        if(data[iWord]>>28 == 1) //TDC header separates events 0x1000XXXX format
                        {
                          unsigned int v1495_header = data[iWord];
                          //printf("TDC header: 0x%x \n", v1495_header);
                          unsigned int t_stop = data[++iWord] & 0xfff;//stop time
                          int v1495_eventID_coda = data[++iWord];//physics event ID recorded from CODA

                          int v1495_eventID_HIGH = data[++iWord];
                          int v1495_eventID_LOW = data[++iWord];

                          //printf("STOP Time = 0x%x => Event ID Coda = 0x%x, from DC HIGH = 0x%x LOW = 0x%x \n",  t_stop, v1495_eventID_coda,v1495_eventID_HIGH, v1495_eventID_LOW );
                          // this doesn't work yet. for codaID = 0x0 it decoes 0x7fff 0xffff for high and low
                          // will skip for now and will use coda event ID for analysis;
                          /*
                          int v1495_eventID = (v1495_eventID_HIGH <<15) + v1495_eventID_LOW; ////physics event ID recorded by Memory card;
                          std::cout << v1495_eventID_HIGH<<15 << "\t "<< v1495_eventID_LOW << "\t" << v1495_eventID << "\n";
                          */

                          //once we got a stop time, we can decode TDC hits:
                          //printf("v1495_board_num = %i \n", v1495_board_num);
                          if(t_stop != 0x2ad){
                            rocs[rocID].tdcs[v1495_board_num].finalizeEvent(codaEventID, v1495_eventID_coda);
                            rocs[rocID].tdcs[v1495_board_num].fillV1495Header(t_stop, 0x0);//0x0 should be replaced with something.
                            rocs[rocID].tdcs[v1495_board_num].fillV1495Hit(v1495_tdc, t_stop);
                          }


/*
                          for(size_t j=0; j<v1495_tdc.size();j++){
                            printf("0x%x \t",v1495_tdc[j]);

                          }
                          printf("\n");
*/
                          v1495_tdc.clear(); //clear tdc hit vector for every new event in the buffer;

                          v1495extraWords = v1495extraWords + 2;
                          i=i+4;

                        }

                        i++;
                      }
                    }


                }
                else if(data[iWord] == 0xe906f018 || data[iWord] == 0xe906f01b) //TW-TDC or QIE
                {
                    unsigned int eventFlag = data[iWord++];
                    if((data[iWord] >> 30) != 0 || (data[iWord] & 0xffff) > 0x0fff)
                    {
                        //ARM dead
                        if(!ARMdead[rocID])
                        {
                            cout << "ARM dead on ROC " << rocID-10 << endl;
                            ARMdead[rocID] = true;
                            ARMdeadFlag = true;
                        }
                        iWord = maxRocWordID;
                        break;
                    }

                    int boardID = ((data[iWord] & 0x0f000000) >> 24) - 9;
                    //cout << "BoardID = " << hex << data[iWord] << "  " << dec << boardID << endl;
                    unsigned int nWordsTDC = data[iWord++] & 0xffff;
                    for(unsigned int i = 0; i < nWordsTDC; ++iWord)
                    {
                        if(data[iWord] == 0xe906e906) continue;
                        if(data[iWord] == nWordsTDC && (i == 0 || i == 1))
                        {
                            ++i;
                        }
                        else if(eventFlag == 0xe906f018)
                        {
                            if((data[iWord] >> 28) == 0) //eventID
                            {
                                rocs[rocID].tdcs[boardID].finalizeEvent(codaEventID, data[iWord]);
                                ++i;
                            }
                            else if((data[iWord] >> 31) != 0) //header
                            {
                                rocs[rocID].tdcs[boardID].fillHeader(data[iWord]);
                                ++i;
                            }
                            else
                            {
                                rocs[rocID].tdcs[boardID].fillHit(data[iWord]);
                                ++i;
                            }
                        }
                        else if(eventFlag == 0xe906f01b)
                        {
                            if((data[iWord] & 0xffff) != 0)
                            {
                                rocs[rocID].tdcs[boardID].finalizeEvent(codaEventID, data[iWord]);
                            }
                            ++i;
                        }
                        else
                        {
                            ++i;
                        }
                    }
                }
                else
                {
                    ++iWord;
                }
            }
        }
        ++codaEventID;
    }

    saveFile->cd();
    saveTree->Write();
    saveFile->Close();

    coda->codaClose();
    return 0;
}
//===========================================================================================

Event::Event()
{
    codaEventID = -1;
    eventID = -1;

    nEntriesExp = 0;
    triggerTime = -1.;

    channels.clear();
    tdcTimes.clear();
}

void Event::setEventID(int cEvtID, int evtID)
{
    codaEventID = cEvtID;
    eventID = evtID;
}

void Event::addHit(unsigned int hit)
{
    int channelID;
    double tdcTime;
    fillInfo(hit, triggerTime, channelID, tdcTime);

    channels.push_back(channelID);
    tdcTimes.push_back(tdcTime);
}

// Added a decoder for V1495 TDC Events=> Ievgen 08/23/2021
void Event::addV1495Hit(std::vector <unsigned int> tdc_word, unsigned int common_stop){

    for(size_t i=0; i< tdc_word.size(); i++){

      int tdc_ch = (int) (tdc_word[i] & 0xff00) >> 8;
      double tdc_time = (double) common_stop - (tdc_word[i] & 0xff);

      channels.push_back(tdc_ch);
      tdcTimes.push_back(tdc_time);
    }
}


void Event::setHeader(unsigned int header)
{
    triggerTime = decodeTime(header);
    nEntriesExp = ((header & 0x0ff00000) >> 20) - 1;
}

void Event::setV1495Header(unsigned int stop_time, unsigned int n_events)
{
    triggerTime = stop_time;
    nEntriesExp = n_events;
}


void Event::fillInfo(unsigned int word, double triggerTime, int& channelID, double& tdcTime)
{
    channelID = ((word & 0xff000000) >> 24) - 0x40;
    tdcTime = triggerTime - decodeTime(word);
    if(tdcTime < 0) tdcTime += 4096.;
}

double Event::decodeTime(unsigned int word)
{
    double fineTime = 4. - (word & 0xf)*4./9.;
    double roughTime = ((word & 0xfff0) >> 4)*4.;

    return roughTime + fineTime;
}

TDC::TDC()
{
    boardID = -1;
    events.clear();
}

void TDC::init()
{
    events.clear();

    Event newEvent;
    events.push_back(newEvent);
}

void TDC::finalizeEvent(int codaEventID, int eventID)
{
    events.back().setEventID(codaEventID, eventID);

    Event newEvent;
    events.push_back(newEvent);
}

void TDC::fillHeader(unsigned int header)
{
    events.back().setHeader(header);
}

void TDC::fillV1495Header(unsigned int stop_time, unsigned int n_events){
    events.back().setV1495Header(stop_time, n_events);
}


void TDC::fillHit(unsigned int hit)
{
    events.back().addHit(hit);
}

void TDC::fillV1495Hit(std::vector <unsigned int> tdc_word, unsigned int common_stop)
{
    events.back().addV1495Hit(tdc_word, common_stop);
}


TString TDC::check()
{
    TString result = "";

    int nErrors1 = 0;  //eventID jump
    int nErrors2 = 0;  //nHits mismatch
    if(events.size() >= 2)
    {
        for(unsigned int i = 1; i < events.size(); ++i)
        {
            //if(events[i].eventID - events[i-1].eventID != 1 && events[i].eventID > 0 && events[i-1].eventID > 0 && events[i].codaEventID - events[i-1].codaEventID < 9000) ++nErrors1;
            if(events[i].eventID - events[i-1].eventID != 1 && events[i].eventID > 0 && events[i-1].eventID > 0 && events[i].codaEventID - events[i-1].codaEventID < 9000) {
               ++nErrors1;
               //cout << i << "  " << events[i].eventID << "  " << events[i-1].eventID << "  " << events[i].codaEventID << "  " << events[i-1].codaEventID << endl;
            }
            if(events[i].tdcTimes.size() != events[i].nEntriesExp && events[i].tdcTimes.size() != 255) ++nErrors2;
            //cout << events[i].tdcTimes.size() << "  " << events[i].nEntriesExp << endl;
        }
    }

    result = Form("%d  %lu  %d  %d", boardID, events.size(), nErrors1, nErrors2);
    return result;
}

void ROC::init()
{
    for(int i = 0; i < nTDCs; ++i) tdcs[i].init();
}

TString ROC::check()
{
    TString result = Form("ROC %02d %d", rocID-10, nTDCs);
    if(tdcs.size() < 1) return result;

    for(int i = 0; i < nTDCs; ++i)
    {
        result = result + " : " + tdcs[i].check();
    }

    return result;
}
