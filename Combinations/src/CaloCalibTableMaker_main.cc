//To run:  g++ -std=c++14 CaloCalibTableMaker.cc -o a.out the ./a.out (add -pthread to use multithreaded)
#include "CaloCalibration/Combinations/inc/CaloCalibTableMaker.hh" 

unsigned int total_ids = 2720 ;
//= 2696 sipms + 16 disk pin diodes + 4 laser pin diodes + 4 spare channels

// inputTable : reads an input archive table
ArchiveTable inputTable(const std::string& archiveFile, const std::string& tag) {
    std::vector<int> roids;
    std::vector<double> peaks;
    std::vector<double> errpeaks;
    std::vector<double> widths;
    std::vector<double> errwidths;
    std::vector<double> sigmas;
    std::vector<double> errsigmas;
    std::vector<double> chisqs;
    std::vector<int> Nhits;

    std::ifstream file(archiveFile);
    if (!file) {
        std::cerr << "[In inputTable()] : Error: Cannot open file " << archiveFile << std::endl;
        exit(1);
    }

    std::string line;
    // Read the file line by line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;

        int roid = 0;
        double peak = 0.;
        double errpeak = 0.;
        double width = 0.;
        double errwidth = 0.;
        double sigma = 0.;
        double errsigma = 0.;
        double chisq = 0.;
        int nhits = 0;

        // Read the first value (roid) separately
        if (std::getline(ss, token, ',')) {
            roid = std::stoi(token);

            // Check the remaining values in the stringstream
            bool all_zeros = true;
            std::string remaining_part;
            std::getline(ss, remaining_part); // Get the rest of the line

            for (char c : remaining_part) {
                if (c != '0' && c != ',' && c != '\n' && c != '\r') {
                    all_zeros = false;
                    break;
                }
            }

            // If all remaining values are '0's, replace them with '1's
            if (all_zeros) {
                std::stringstream modified_ss(remaining_part);
                std::string zero_token;
                // Discard the first value, then process the rest
                while (std::getline(modified_ss, zero_token, ',')) {
                     // Reset all 0,0,0,0 etc. to 1s
                    peak = 1.0;
                    errpeak = 1.0;
                    width = 1.0;
                    errwidth = 1.0;
                    sigma = 1.0;
                    errsigma = 1.0;
                    chisq = 1.0;
                    nhits = 1;
                    
                    peaks.push_back(peak);
                    errpeaks.push_back(errpeak);
                    widths.push_back(width);
                    errwidths.push_back(errwidth);
                    chisqs.push_back(chisq);
                    sigmas.push_back(sigma);
                    errsigmas.push_back(errsigma);
                    Nhits.push_back(nhits);
                }
                
                roids.push_back(roid);
            } else {
                ss.clear();
                ss.str(remaining_part);

                // Extract each token separated by a comma
                if (std::getline(ss, token, ',')) peak = std::stod(token);
                if (std::getline(ss, token, ',')) errpeak = std::stod(token);
                if (std::getline(ss, token, ',')) width = std::stod(token);
                if (std::getline(ss, token, ',')) errwidth = std::stod(token);
                if (std::getline(ss, token, ',')) sigma = std::stod(token);
                if (std::getline(ss, token, ',')) errsigma = std::stod(token);
                if (std::getline(ss, token, ',')) chisq = std::stod(token);
                if (std::getline(ss, token, ',')) nhits = std::stoi(token);

                // Populate the vectors
                roids.push_back(roid);
                peaks.push_back(peak);
                errpeaks.push_back(errpeak);
                widths.push_back(width);
                errwidths.push_back(errwidth);
                chisqs.push_back(chisq);
                sigmas.push_back(sigma);
                errsigmas.push_back(errsigma);
                Nhits.push_back(nhits);
            }
        }
    }

    return ArchiveTable(roids, peaks, errpeaks, widths, errwidths, sigmas, errsigmas, chisqs, Nhits, tag);
}
// combineAlg : called by main function, this funciton does the hard work
RecoTable combineAlg(std::vector<ArchiveTable> tables){
  std::vector<int> roids;
  std::vector<double> peaks;
  // Fill ROID:
  for(unsigned int i = 0; i < total_ids; i++){
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
    if(id > 2695) peaks.push_back(1);
    else peaks.push_back(average_peak/tables.size());
  }
   
  RecoTable outputs(roids, peaks);
  return outputs;
}


void WriteOutput(const RecoTable& table) {
  std::ofstream fw("RecoTable.txt", std::ofstream::out);

  // Check if the file was successfully opened
  if (!fw) {
    std::cerr << "Error: Unable to open RecoTable.txt for writing." << std::endl;
    return;
  }

  fw << "TABLE CalEnergyCalib" << "\n";
  
  if (!table.roid_.empty()) {
      for (size_t i = 0; i < table.roid_.size(); ++i) {
          // Use '\n' instead of std::endl for better performance
          fw << table.roid_[i] << "," << table.ADC2MeV_[i] << "\n";
      }
  }

  // Check for write errors after all operations are complete
  if (!fw) {
    std::cerr << "Error: An error occurred while writing to RecoTable.txt." << std::endl;
  }
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
