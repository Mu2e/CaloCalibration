//To run:  g++ -std=c++14 CaloCalibTableMaker.cc -o a.out the ./a.out (add -pthread to use multithreaded)
#include "CaloCalibration/Combinations/inc/CaloCalibTableMaker.hh" 


// inputTable : reads an input archive table
ArchiveTable inputTable(const char *ArchiveFile, const char *tag){
  std::vector<int> roids;
  std::vector<double> peaks;
  std::vector<double> errpeaks;
  std::vector<double> widths;
  std::vector<double> errwidths;
  std::vector<double> sigmas;
  std::vector<double> errsigmas;
  std::vector<double> chisqs;
  std::vector<int> Nhits;
  
  FILE *fC = fopen(ArchiveFile, "r");
  
	if ( fC == NULL) {
		std::cout<<"[In inputTable()] : Error: Cannot open file "<<std::endl;
		exit(1);
	}
	int roid, nhits;
	float peak, errpeak, width, errwidth, sigma, errsigma, chisq;
	
  while(fscanf(fC, "%i,%f,%f,%f,%f,%f,%f,%f,%i\n", &roid, &peak, &errpeak, &width, &errwidth, &sigma, &errsigma, &chisq, &nhits)!=EOF){
    roids.push_back(roid);
	  peaks.push_back(peak);
	  errpeaks.push_back(errpeak);
	  widths.push_back(width);
	  errwidths.push_back(errwidth);
	  chisqs.push_back(chisq);
	  sigmas.push_back(sigma);
	  errsigmas.push_back(errsigma);
	  Nhits.push_back(nhits);
	  std::cout<<roid<<","<<peak<<","<<errpeak<<","<<width<<","<<errwidth<<","<<sigma<<","<<errsigma<<","<<chisq<<","<<nhits<<std::endl;
  }
  
	fclose(fC);
	ArchiveTable input(roids, peaks, errpeaks, widths, errwidths, sigmas, errsigmas, chisqs, Nhits, tag);
  return input;
}

// combineAlg : called by main function, this funciton does the hard work
RecoTable combineAlg(std::vector<ArchiveTable> tables){
  std::vector<int> roids;
  std::vector<double> peaks;
  // Fill ROID:
  for(unsigned int i = 0; i < 1348*2; i++){ // FIXME - needs to use CaloID class and also no hardcoding
    roids.push_back(i);
  }
  // Loop over ROID, for each table find ID value and average, fill new table
  for(auto const& id: roids){
    double average_peak = 0;
    for(unsigned int i = 0; i < tables.size() ; i++){
      
      if(tables[i].tag_ == "mip"){
        if(tables[i].peak_[id] != 0) average_peak += 20/tables[i].peak_[id];
        else average_peak += 0;
      } else {
        average_peak += tables[i].peak_[id]; //assumes source is already ADC2MeV
      }
    }
    peaks.push_back(average_peak/tables.size());
  }
   
  RecoTable outputs(roids, peaks);
  return outputs;
}


void WriteOutput(RecoTable table) {
  ofstream fw("RecoTable.txt", std::ofstream::out);
  if (fw.is_open())
  {
    fw <<"TABLE CalEnergyCalib"<<"\n";
    for(unsigned int i = 0; i< table.roid_.size(); i++){
      fw << table.roid_[i]<<","<<table.ADC2MeV_[i]<<"\n";
    }
  }
  fw.close();
}

void removeFirstLine(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    std::ofstream tempFile("temp.txt"); // Create a temporary file
    if (!tempFile.is_open()) {
        std::cerr << "Error: Could not create temporary file." << std::endl;
        inputFile.close();
        return;
    }
    std::string line;
    // Read and discard the first line
    std::getline(inputFile, line); 

    // Read and write the rest of the lines
    while (std::getline(inputFile, line)) {
        tempFile << line << std::endl;
    }
    inputFile.close();
    tempFile.close();
    // Delete the original file
    if (std::remove(filename.c_str()) != 0) {
        std::cerr << "Error: Could not delete original file." << std::endl;
        return;
    }
    // Rename the temporary file to the original file's name
    if (std::rename("temp.txt", filename.c_str()) != 0) {
        std::cerr << "Error: Could not rename temporary file." << std::endl;
        return;
    }
    std::cout << "First line removed successfully from " << filename << std::endl;
}

int main(int argc, char ** argv){
  bool source = argv[2];
  bool MIP = argv[1];
  std::vector<ArchiveTable> tables;
  // get source table
  if(source){
    removeFirstLine("Source.txt");
    ArchiveTable sourceInput = inputTable("Source.txt", "source");
    tables.push_back(sourceInput);
  }
  // get MIP table
  if(MIP){
    removeFirstLine("Cosmics.txt");
    ArchiveTable mipInput = inputTable("Cosmics.txt", "mip");
    tables.push_back(mipInput);
  }
  // pass to combination:
  RecoTable output = combineAlg(tables);
  //write out CSV:
  WriteOutput(output);
  return 0;
}
