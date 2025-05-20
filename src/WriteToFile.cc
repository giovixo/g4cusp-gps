#include "WriteToFile.hh"

WriteToFile::WriteToFile() {
    outputFile.open("output.txt");
}

WriteToFile::~WriteToFile() {
    outputFile.close();
}

void WriteToFile::print(const std::string& content) {
    outputFile << content << std::endl;
}