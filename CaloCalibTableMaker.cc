//To run:  g++ -std=c++14 CaloCalibTableMaker.cc -o a.out the ./a.out (add -pthread to use multithreaded)
#include "CaloCalibTableMaker.hh" 


// inputTable : reads an input archive table
Table inputTable(const char *ArchiveFile, const char *tag){
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
	Table input(roids, peaks, errpeaks, widths, errwidths, chisqs, tag);
  return input;
}
// outputTable : returns final table in .txt form
void outputTable(Table table){

}

// combineAlg : called by main function, this funciton does the hard work
Table combineAlg(std::vector<Table> tables){
  Table output;
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
  
  for(unsigned int i = 0; i < tables.size() ; i++){
    std::cout<<"looking at table "<<tables[i].algtag_<<std::endl;
    roids.push_back(tables[i].roid_);
    
    
  }
  Table output(roids, peaks, errpeaks, widths, errwidths, chisqs, "reco combined");
}

int main(){
  std::vector<Table> tables;
  // get source table
  Table sourceInput = inputTable("testdata.csv", "source");
  tables.push_back(sourceInput);
  // get MIP table
  Table mipInput = inputTable("testdata.csv", "mip");
  tables.push_back(mipInput);
  // get laser table
  Table laserInput = inputTable("testdata.csv", "laser");
  tables.push_back(laserInput);
  // pass to combination:
  combineAlg(tables);
  return 0;
}
