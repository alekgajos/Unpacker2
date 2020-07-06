#include <iostream>
#include <map>
#include <vector>
#include <string>

#include <TTree.h>
#include <TFile.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Unpacker2D.h"
#include "EventIII.h"
#include "TDCChannel.h"


//############
//
// DEFINE THE NUMBER OF ENDPOINTS
//
//############

#define ENDPOINTS 2


using namespace std;

UInt_t ReverseHexDJ(UInt_t n) {
	UInt_t a, b, c, d, e;
	a = n & 0x000000ff;
	b = n & 0x0000ff00;
	c = n & 0x00ff0000;
	d = n & 0xff000000;

	a <<= 8;
	b >>= 8;
	c <<= 8;
	d >>= 8;

	e = a|b|c|d;

	return e;
}

UInt_t ReverseHexTDC(UInt_t n) {
	UInt_t a, b, c, d, e;
	a = n & 0x000000ff;
	b = n & 0x0000ff00;
	c = n & 0x00ff0000;
	d = n & 0xff000000;

	a <<= 24;
	b <<= 8;
	c >>= 8;
	d >>= 24;

	e = a|b|c|d;

	return e;
}

void Unpacker2D::BuildEvent(EventIII* e, map<UInt_t, vector<UInt_t> >* m, map<UInt_t, double>* refTimes) {
	UInt_t data;
	UInt_t fine;
	UInt_t coarse;
	UInt_t rising;
	double refTime = 0;
	double time = 0;
	double time_t = 0;

	map<UInt_t, vector<UInt_t> >::iterator m_it;
	for (m_it = m->begin(); m_it != m->end(); m_it++) {
		TDCChannel* tc = e->AddTDCChannel(m_it->first);
		for (int i = 0; i < m_it->second.size(); i++) {

			data = m_it->second[i];

			fine = data & 0xff;
			coarse = (data >> 8) & 0xffff;
			rising = (data >> 31);

			if (useTDCcorrection == true)
				time = (coarse * 2.5) - ((TDCcorrections[m_it->first]->GetBinContent(fine + 1)) / 1000.0);
			else
				time = (coarse * 2.5) - (fine * 0.0208333333);
				
			refTime = refTimes->find((int)((m_it->first - 2100) / 105) * 105 + 2100)->second;

			if (rising == 0) {
				if (time - refTime < 0)
					time_t = time + ((0x10000 * 2.5) - refTime);
				else 
					time_t = time - refTime;

				tc->AddLead(time_t);
			}
			else {
				if (time - refTime < 0)
					time_t = time + ((0x10000 * 2.5) - refTime);
				else
					time_t = time - refTime;

				tc->AddTrail(time_t);
			}
	
		}
	}
}

void Unpacker2D::ParseConfigFile(string f, string s) {
	  // parsing xml config file
  boost::property_tree::ptree tree;

  try {
    boost::property_tree::read_xml(s, tree);
  } catch(boost::property_tree::xml_parser_error e) {
    cerr << "ERROR: Failed to read config file" << endl;
    exit(0);
  }

  // get the config options from the config file
  try {
    if (tree.get<string>("READOUT.DEBUG") == "ON") {
      debugMode = true;
    }
  } catch (exception e) {
    cerr << "ERROR: Incorrect config file structure" << endl;
    exit(0);
  }

  if (debugMode == true)
    cerr<<"DEBUG mode on"<<endl;

  // get the first data source entry in the config file
  boost::property_tree::ptree readoutTree = tree.get_child("READOUT");
  string type;
  string address_s;
  UInt_t address;
  string hubAddress;
  string correctionFile;
  int channels = 0;
  int offset = 0;
  //int resolution = 0;
  //int referenceChannel = 0;
  string measurementType("");
  highest_channel_number = 0;

  // iterate through entries and create appropriate unpackers
  for (const auto& readoutEntry : readoutTree) {
    // read out values from xml entry
    if ((readoutEntry.first) == "DATA_SOURCE") {
      type = (readoutEntry.second).get<string>("TYPE");
      address_s = (readoutEntry.second).get<string>("TRBNET_ADDRESS");
      hubAddress = (readoutEntry.second).get<string>("HUB_ADDRESS");
      //referenceChannel = (readoutEntry.second).get<int>("REFERENCE_CHANNEL");
      correctionFile = (readoutEntry.second).get<string>("CORRECTION_FILE");

      if (correctionFile.compare("raw") != 0){
        cerr << "WARNING: The TDC correction file path was set in the XML config file of the Unpacker!" << endl;
        cerr << "This file path should be defined in the user parameters JSON file instead." << endl;
        cerr << "The setting from the XML file fill be ignored!" << endl;
      }

      if (type == "TRB3_S" || type == "DJPET_ENDP") {

        // create additional unpackers for internal modules
        boost::property_tree::ptree modulesTree = (readoutEntry.second).get_child("MODULES");
        for (const auto& module : modulesTree) {
          type = (module.second).get<string>("TYPE");
          address_s = (module.second).get<string>("TRBNET_ADDRESS");
          address = std::stoul(address_s, 0 , 16);
          channels = (module.second).get<int>("NUMBER_OF_CHANNELS");
          offset = (module.second).get<int>("CHANNEL_OFFSET");
          //resolution = (module.second).get<int>("RESOLUTION");
          measurementType = (module.second).get<string>("MEASUREMENT_TYPE");

          tdc_offsets[address] = offset;
          if( offset + channels > highest_channel_number ){
            highest_channel_number = offset + channels;
          }
        }
      } else {
        cerr << "Incorrect configuration in the xml file!" << endl;
        cerr << "The DATA_SOURCE entry is missing!" << endl;
      }
    }
  }

}

void Unpacker2D::UnpackSingleStep(const char* hldFile, const char* configFile, int numberOfEvents, int refChannelOffset, const char* TOTcalibFile, const char* TDCcalibFile) {
	eventsToAnalyze = numberOfEvents;
	this->refChannelOffset = refChannelOffset;

	//TODO: read config file here
	//tdc_offsets[0xa110] = 0;
	ParseConfigFile(string(hldFile), string(configFile));

	useTDCcorrection = false;
	if (strlen(TDCcalibFile) != 0) {
		useTDCcorrection = loadTDCcalibFile(TDCcalibFile);
		printf("CAL LOADED: %d\n", useTDCcorrection);
	}

	DistributeEventsSingleStep(string(hldFile));
}

bool Unpacker2D::loadTDCcalibFile(const char* calibFile) {
	TFile* f = new TFile(calibFile, "READ");
	TDirectory* dir = gDirectory->GetDirectory("Rint:/");
	
	if (f->IsOpen()) {
		TDCcorrections = new TH1F*[highest_channel_number];
		for (int i=2100;i<highest_channel_number;i++) {
			TH1F* tmp = dynamic_cast<TH1F*>(f->Get(Form("correction%d", i)));

			if (tmp) {
				TDCcorrections[i] = dynamic_cast<TH1F*>(tmp->Clone(tmp->GetName()));
				TDCcorrections[i]->SetDirectory(dir);
			} else {
				TDCcorrections[i] = nullptr;
			}
		}

		if(debugMode){
		      cerr << "Loaded TDC nonlinearity corrections." << endl;
    		}

  	}else{
    		if(debugMode){
      			cerr << "The TDC calibration file " << calibFile << " could not be properly opened." << endl;
			cerr << "TDC nonlinearity correction will not be used!" << endl;
    		}	

    		f->Close();
		delete f;

		return false;
	}
	f->Close();
	delete f;
	return true;
}

void Unpacker2D::DistributeEventsSingleStep(string filename) {
	ifstream* file = new ifstream(filename.c_str());
	

	if (file->is_open()) {
	
		EventIII* eventIII = new EventIII();
		
		string newFileName = filename + ".root";
		TFile* newFile = new TFile(newFileName.c_str(), "RECREATE");
		TTree* newTree = new TTree("T", "Tree");

		newTree->Branch("eventIII", "EventIII", &eventIII, 64000, 99);

		TH1F* h_errors = new TH1F("h_errors", "h_errors", 10, 0, 10);
		//TH1F* h_ts_width = new TH1F("h_ts_width", "h_ts_width", 500000, 45000, 55000);
		TH1F* h_ts_width = new TH1F("h_ts_width", "h_ts_width", 10000, 0, 1000000);
		TH1F* h_ts_diff = new TH1F("h_ts_diff", "h_ts_diff", 100000, -250, 250);
		TH1F* h_coarse = new TH1F("h_coarse", "h_coarse", 65535, 0, 65535);
		TH1F* h_fine = new TH1F("h_fine", "h_fine", 128, 0, 128);

		UInt_t data4;
		Int_t nBytes;
		UInt_t nEvents = 0;

		UInt_t queueSize = 0;
		UInt_t queueDecoding = 0;
		UInt_t subSize = 0;
		UInt_t subDecoding = 0;
		UInt_t ftabId = 0;
		UInt_t ftabSize = 0;
		UInt_t ftabTrgn = 0;
		UInt_t ftabDbg = 0;
		UInt_t dataCtr = 0;
		UInt_t channel = 0;
		UInt_t currentOffset = 0;

		UInt_t ftabWords = 0;

		UInt_t badIdCtr = 0;

		map<UInt_t, UInt_t>::iterator offsets_it;
		map<UInt_t, vector<UInt_t> > tdc_channels;
		map<UInt_t, double> refTimes;
		map<UInt_t, double> refTimes_previous;

		UInt_t prevTrgId = 0;

		bool trgSequence = false;
		bool trgSequenceHold = false;

		int built = 0;
		int skip = 0;

		bool missing_ref = false;
		int refTimesCtr = 0;
		int refTimesCtrPrevious = 0;

		int missingRefCtr = 0;

		// skip the first entry
		file->ignore(32);

		unsigned int packet_buf[1024];
		int packet_buf_ctr = 0;

		while(!file->eof()) {
			nBytes = 0;
			if (nEvents % 10000 == 0) {
				printf("%d\n", nEvents);
			}
// if (nEvents == 5140)
// 	break;
			// printf("\nEvent no:%d\n", nEvents);

			eventIII->Clear();

			// queue headers
			// size
			file->read((char*) & data4, 4);
			nBytes += 4;
			queueSize = data4 / 4;

			nEvents++;
			

			// printf("queue size: %d x:%08x\n", queueSize, data4);

			// skip bad entries
			if (queueSize < 0x10) {
				file->ignore(28);
				nBytes += 28;
				// if (debugMode)
					// printf("Skipping too small queue\n");
					h_errors->Fill(2);
					newTree->Fill();
				continue;
			}

//			printf("processing\n");

			// decoding
			file->read((char*)&data4, 4);
			nBytes += 4;
			queueDecoding = data4;

			// printf("%08x\n",data4);

			// skip some headers
			file->ignore(8);
			nBytes += 8;

			// subevent
			// sub size
			file->read((char*)&data4, 4);
			nBytes += 4;
			subSize = data4 & 0xffff;

			// printf("%08x\n",data4);

			// sub decoding
			file->read((char*)&data4, 4);
			nBytes += 4;
			subDecoding = data4;

			// printf("%08x\n",data4);

			// skip some headers
			file->ignore(24); 
			nBytes += 24;

			queueSize -= 12;

			//printf("qs before endp %d\n", queueSize);

			// refTimes.clear();

			refTimesCtrPrevious = refTimesCtr;
			refTimesCtr = 0;

			packet_buf_ctr = 0;

			prevTrgId = ftabTrgn;
			
			while(!file->eof()) {
				// ftab header
				// ftab size and id
				file->read((char*)&data4, 4);
				nBytes += 4;
				
				// printf("size and id: %x\n", data4);
				
				data4 = ReverseHexDJ(data4);
				ftabSize = data4 >> 16;
				ftabId = data4 & 0xffff;

				ftabWords = ftabSize - 2;

				// printf("%08x\n",data4);

				// ftab trigger number and debug
				file->read((char*)&data4, 4);
				nBytes += 4;
				data4 = ReverseHexDJ(data4);
				ftabDbg = data4 >> 16;
				ftabTrgn = data4 & 0xffff;

				queueSize -= 2;

				// printf("%08x\n",data4);
				// printf("id: %x trg: %d size: %d, fWords: %d\n", ftabId, ftabTrgn, ftabSize, ftabWords);


				offsets_it = tdc_offsets.find(ftabId); 
				if (offsets_it == tdc_offsets.end()) {
					// printf("Wrong ftab ID: %04x, evt no: %d\n", ftabId, nEvents);
					h_errors->Fill(3);
					file->ignore( (ftabSize - 2) * 4);
					nBytes += (ftabSize - 2) * 4;
					queueSize -= ftabSize - 2;
					// printf("%d %d %d\n", nBytes, queueSize, ftabSize);
					badIdCtr++;
					trgSequence = false;
					if (queueSize == 0) { break; }
					else if (queueSize == 1) { file->ignore(4); nBytes += 4; queueSize -= 1; break; }
					continue;
				}
				currentOffset = offsets_it->second;
				// std::cerr<<std::hex<<ftabId<<std::dec<<" "<<ftabSize<<" "<<currentOffset<<std::endl;
				if (nEvents > 0) refTimes_previous[currentOffset] = refTimes[currentOffset];

				// ftab data
				while(!file->eof()) {
					file->read((char*)&data4, 4);
					nBytes += 4;
					queueSize--;
					data4 = ReverseHexTDC(data4);

					ftabWords--;

					// printf("%08x  fw: %d\n",data4, ftabWords);

					// if (data4 == 0xffffffff) {
					// 	file->read((char*)&data4, 4);
					// 	printf("last word: %x, q: %d, fs: %d, fw: %d\n", data4, queueSize, ftabSize, ftabWords);
					// 	nBytes += 4;
					// 	queueSize--;
					// 	break;
					// }

					if (ftabWords == 1) { 
						file->read((char*)&data4, 4);
						// printf("last word: %x, q: %d, fs: %d, fw: %d\n", data4, queueSize, ftabSize, ftabWords);
						nBytes += 4;
						queueSize--;
						break; 
					}

					if ((data4 >> 24) != 0xfc) {
						channel = (data4 >> 24) & 0x7f;

						h_coarse->Fill((data4 >> 8) & 0xffff);
						h_fine->Fill(data4 & 0xff);

// if (channel != 104)
						// printf("%d %d\n", channel, (data4) & 0x80);
					
						//tdc_channels[channel].push_back(data4);

//printf("%x\n", data4);
						
						if (channel == 104) {
							if (useTDCcorrection == true) {
								//if (nEvents > 0) refTimes_previous[currentOffset] = refTimes[currentOffset];
								refTimes[currentOffset] = (((data4 >> 8) & 0xffff) * 2.5) - ((TDCcorrections[channel + currentOffset]->GetBinContent((data4 & 0xff) + 1)) / 1000.0);

								refTimesCtr++;

							}
							else {
								//if (nEvents > 0) refTimes_previous[currentOffset] = refTimes[currentOffset];
								refTimes[currentOffset] = (((data4 >> 8) & 0xffff) * 2.5) - ((data4 & 0xff) * 0.0208333333);

								refTimesCtr++;
							}

						}
						else {
							// if (channel != 99)
								// printf("%08x\n", data4);
								tdc_channels[channel + currentOffset].push_back(data4);
						}

						//tdc_channels[channel].push_back(data4);
					}

					dataCtr++;

				}

				if (queueSize < 4) {
					if (queueSize == 1)
						file->ignore(4);

					// printf("refTimes %d %d\n", refTimes.size(), missing_ref, refTimesCtr);

		// printf("No: %d trg diff: %d  ftrgnr: %d sequence: %d hold: %d refsctr: %d, ref: %lf", nEvents, (ftabTrgn - prevTrgId), ftabTrgn, trgSequence, trgSequenceHold, refTimesCtr, refTimes_previous[2100]);
					if (refTimesCtr == ENDPOINTS) {
					// if (refTimesCtr == 4) {
					//if (refTimesCtr == 1) {
						// printf("No: %d trg diff: %d   ref: %lf\n", nEvents, (ftabTrgn - prevTrgId), refTimes_previous[2100]);
						if (nEvents > 0) {
							if (missing_ref == false && (ftabTrgn - prevTrgId) == 1 && trgSequence == true && trgSequenceHold == false) {
								h_ts_diff->Fill((refTimes[2100] - refTimes[2205]) - (refTimes_previous[2100] - refTimes_previous[2205]));
								h_ts_width->Fill(refTimes[2100] - refTimes_previous[2100]);
								// h_ts_width->Fill(refTimes[2205] - refTimes_previous[2205]);

								map<UInt_t, double>::iterator ref_it = refTimes_previous.begin();
								while( ref_it != refTimes_previous.end()) {
									TDCChannel* tc = eventIII->AddTDCChannel(channel + ref_it->first);
									tc->AddLead(refTimes_previous[ref_it->first]);
									// TDCChannel* tc2 = eventIII->AddTDCChannel(channel + 2205);
									// tc2->AddLead(refTimes_previous[2205]);
									ref_it++;
								}

								// printf("%d - %d: %lf %lf\n", nEvents, ftabTrgn, refTimes[2100], abs(refTimes[2100] - refTimes_previous[2100]));
// printf("\tbuilding\n");
built++;
								// if ( abs(refTimes[2100] - refTimes_previous[2100]) < 50100.0 || abs(refTimes[2100] - refTimes_previous[2100]) > 40900.0) {
								BuildEvent(eventIII, &tdc_channels, &refTimes_previous);
								// }
								// newTree->Fill();
								h_errors->Fill(0); // event ok
								// refTimes.clear();

								trgSequence = true;
							}
							else {
								skip++;
								// printf("\tskipping\n");
							}							
							// newTree->Fill();
							missing_ref = false;
						}

						if ((ftabTrgn - prevTrgId) == 1) {
							trgSequence = true;
						}
						else if((ftabTrgn - prevTrgId) != 1) {
							trgSequence = false;
						}

						if (refTimesCtrPrevious == ENDPOINTS) {
							trgSequenceHold = false;
						}
					}
					else  {
						skip++;
						// printf("\tskipping\n");
						// if (debugMode)
							// printf("Missing reference time\n");
						// newTree->Fill();
						h_errors->Fill(1); // missing ref time
						missing_ref = true;

						missingRefCtr++;

						trgSequenceHold = true;
						trgSequence = false;
					}

					

					tdc_channels.clear();


					break;
				}

			}

			newTree->Fill();
			// nEvents++;

			if (nEvents == eventsToAnalyze) {
				printf("Max timeslots reached\n");
				break;
			}
		}
		h_fine->Write();
		h_coarse->Write();
		h_errors->Write();
		newFile->Write();
		delete newTree;
		file->close();

		// printf("bad id ctr: %d\n", badIdCtr);
		// printf("built: %d skip: %d\n", built, skip);
		// printf("missing ref ctr: %d\n", missingRefCtr);
	}
	else { cerr<<"ERROR: failed to open data file"<<endl; }
}
