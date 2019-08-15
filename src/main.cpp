#include <TROOT.h>
#include <string>
#include <iostream>
#include <TApplication.h>
#include <chrono>
#include "gwmSimulation.h"

using namespace std;

int main(int argc, char* argv[]) {
  try {
    if(argc != 2) {
      string err = "Incorrect number of arguments! Needs name of output file";
      throw(err);
    } else {
      simulation SIM;
      auto start = chrono::high_resolution_clock::now();
      SIM.run(argv[1]);
      auto stop = chrono::high_resolution_clock::now();
      auto duration = chrono::duration_cast<chrono::seconds>(stop-start);
      cout<<"Elapsed time: "<<duration.count()<<" seconds"<<endl;
    }
  } catch(string err) {
    cout<<err<<endl; 
  }
}
