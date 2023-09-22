//To run:  g++ -std=c++14 CaloCalibTableMaker.cc -o a.out the ./a.out (add -pthread to use multithreaded)
#include "CaloCalibration/SourceCalib/inc/CaloCalibTableMaker.hh" 


// inputTable : reads an input archive table
ArchiveTable inputTable(const char *ArchiveFile, const char *tag){
  std::vector<int> roids;
  std::vector<double> peaks;
  std::vector<double> widths;
  std::vector<double> errpeaks;
  std::vector<double> errwidths;
  std::vector<double> chisqs;
  
  FILE *fC = fopen(ArchiveFile, "r");
	if ( fC == NULL) {
		std::cout<<"[In inputTable()] : Error: Cannot open file "<<std::endl;
		exit(1);
	}
	int roid;
	float peak, errpeak, width, errwidth, chisq;
	
  while(fscanf(fC, "%i,%f,%f,%f,%f,%f\n", &roid, &peak, &errpeak, &width, &errwidth, &chisq)!=EOF){
    roids.push_back(roid);
	  peaks.push_back(peak);
	  errpeaks.push_back(errpeak);
	  widths.push_back(width);
	  errwidths.push_back(errwidth);
	  peaks.push_back(peak);
	  chisqs.push_back(chisq);
  }
  
	fclose(fC);
	ArchiveTable input(roids, peaks, errpeaks, widths, errwidths, chisqs, tag);
  return input;
}

// combineAlg : called by main function, this funciton does the hard work
RecoTable combineAlg(std::vector<ArchiveTable> tables){
  std::cout<<" combining tables "<<std::endl;
  ArchiveTable output;
  std::vector<int> roids;
  std::vector<double> peaks;
  std::vector<double> widths;
  std::vector<double> errpeaks;
  std::vector<double> errwidths;
  std::vector<double> chisqs;
  double p_av = 0;
  double errp_av = 0;
  double w_av = 0;
  double errw_av = 0;
  double chisq_av = 0; //TODO - how to deal with this?
  // Fill ROID:
  for(unsigned int i = 0; i < 1348*2; i++){ // FIXME - needs to use CaloID class and also no hardcoding
    roids.push_back(i);
  }
  // Loop over ROID, for each table find ID value and average, fill new table
  for(auto const& id: roids){
    double average_peak = 0;
    for(unsigned int i = 0; i < tables.size() ; i++){
      //std::cout<<"looking at table "<<tables[i].algtag_<<std::endl;
      average_peak += tables[i].peak_[id];
    }
    peaks.push_back(average_peak/tables.size());
    errpeaks.push_back(0);//TODO
  }
   
  RecoTable outputs(roids, peaks, errpeaks, 1);
  return outputs;
}


void WriteOutput(RecoTable table) {
  ofstream fw("RecoTable.txt", std::ofstream::out);
  if (fw.is_open())
  {
    for(unsigned int i = 0; i< table.roid_.size(); i++){
      fw << table.roid_[i]<<","<<table.peak_[i]<<","<<table.errpeak_[i]<<"\n";
    }
  }
  fw.close();
}

int main(int argc, char ** argv){
  bool source = argv[1];
  bool MIP = argv[2];
  bool laser = argv[3];
  std::vector<ArchiveTable> tables;
  // get source table
  if(source){
    ArchiveTable sourceInput = inputTable("CaloCalibration/data/sourcedata_test.csv", "source");
    tables.push_back(sourceInput);
  }
  // get MIP table
  if(MIP){
    ArchiveTable mipInput = inputTable("CaloCalibration/data/cosmicdata_test.csv", "mip");
    tables.push_back(mipInput);
  }
  // get laser table
  if(laser){
    ArchiveTable laserInput = inputTable("CaloCalibration/data/laserdata_test.csv", "laser");
    tables.push_back(laserInput);
  }
  // pass to combination:
  RecoTable output = combineAlg(tables);
  //write out CSV:
  WriteOutput(output);
  return 0;
}
