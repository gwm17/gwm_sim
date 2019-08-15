/*SRIMscript.cpp
 *Takes a SRIM file and cleans it for use with anasen analysis. Mostly just makes everything have uniform units,
 *removes headers and unecessary info
 *
 *Maria A. -- a long time ago 
 *
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>

using namespace std;

int main(int argc, char **argv)
{
    string alouminum(argv[1]);//original filename
    
    vector<string> s2;//will store keV/MeV (unneccessary)
    vector<double> ionE,dE_dx_e, dE_dx_n, Range, longStrag, latStrag;//will store values to save.
    string keV="keV";//this is just here to utilize in comparisons;
    
    ifstream in(alouminum.c_str());//open in file
    string input/*temporary*/,Straggling="Straggling"/*for comparison*/,commentline/*to save comment line*/;
    //This while loop reads the crap before the actual data.
    while (in>>input)//This saves in input automatically. if EOF will kill it (will never EOF)
    {
        if (input==Straggling)//Check if "Straggling is reached.
        {
            string junk;
            in>>input;//Second Straggling
            in>>commentline;//save this for later.
            getline(in, junk); //toss line of ---'s
            break;//kill loop since we can start reading data.
        }
        input.clear();
    }
    //Reads the rest until it reaches  a line that begins with '-'
    double ie;
    float dedxe;
    float dedxn;
    float range;
    float longstrg;
    float latstrg;
    string units;
    string junk;
    while(in>>ie)
    {
        in>>units;
        if(units=="keV") ie = ie/1000.0;
        
        in>>dedxe;
        in>>dedxn;
        ionE.push_back(ie);
        dE_dx_e.push_back(dedxe);
        dE_dx_n.push_back(dedxn);
        
        in>>range;
        in>>units;
        if(units=="mm") {range = range/10.0;}
        else if (units=="um") {range=range/10000.0;}
        else if (units=="m") {range=range*100;}
        Range.push_back(range);
 
        in>>longstrg;
        in>>units;
        if(units=="mm") {longstrg = longstrg/10.0;}
        else if (units=="um") {longstrg=longstrg/10000.0;}
        else if (units=="m") {longstrg=longstrg*100;}
        longStrag.push_back(longstrg);
 
        in>>latstrg;
        in>>units;
        if(units=="mm") {latstrg = latstrg/10.0;}
        else if (units=="um") {latstrg=latstrg/10000.0;}
        else if (units=="m") {latstrg=latstrg*100;}
        latStrag.push_back(latstrg);
    }
     
    string outputfilename=alouminum.substr(0, alouminum.find("."));//output file name
    outputfilename += ".eloss";
    ofstream out(outputfilename.c_str());//open file.
    
    for(unsigned int i=0; i<ionE.size();i++)//save data.
    {
        out<<fixed<<showpoint;
        out<<setprecision(8);
        out<<setw(10)<<ionE[i]<<"\t"<<dE_dx_e[i]<<"\t"<<dE_dx_n[i]<<"\t"<<Range[i]<<"\t"
           <<longStrag[i]<<"\t"<<latStrag[i]<<endl;
    }
    
	return 0;
}
