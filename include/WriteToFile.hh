#ifndef WriteToFile_H
#define WriteToFile_H 1

#include <fstream>

class  WriteToFile {
public:
    WriteToFile();
    ~WriteToFile();

public:
    void print(const std::string& content);

private:
    std::ofstream outputFile;
};

#endif