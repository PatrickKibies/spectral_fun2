#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <fftw3.h>


/*
 * Small example for using the FFTW3 library in order to perform a 1D discrete fourier transform from
 * real number input to complex output. Including evaluaton of phase for output.
 * Author: Patrick Kibies
 * Date: 2014-10-01
 */


using namespace std;


void printCoeffs ( int i, double c, double s ) //print Fourier coefficients
{
    double a = sqrt(c*c + s*s); //amplitude by pythagoras
    double phi = 0; //phase angle
    
    printf ( "%5d\t%13.6f %13.6f %13.6f", i, c, s, a );
    //if (a>0.000001) { //for non negligible amplitudes
	phi = atan2(s,c) * ( 180 / M_PI ); //calculate phase in degree
        printf ( " %10.3f", phi );
    //}
    printf("\n");
}

void printFrequencies ( int n, fftw_complex *out ) //n is number of samples
{
    int outcount = n; //outcount is number of fourier coefficients
    outcount /=2;
    outcount +=1;
  
  
    int i = 0;
    printf("%5s\t%13s %13s %13s %10s\n","i", "cos", "sin", "Ampl.", "Phase" ); //output table header
    
    printCoeffs (0, out[0][0]/n, 0 ); 
    for ( i=1; i<outcount; i++ ) 
	printCoeffs (i, 2*out[i][0]/n, -2*out[i][1]/n);
    if ( n % 2 == 0 )
	printCoeffs (i, out[i][0]/n, 0 );
}

int main(int argc, char **argv) {
  unsigned int sample_count = 0; //number of samples
  string STRING, delimiter, token; //who needs strings?
  vector<string> lines; //input holder

  fftw_plan p; //no plan? change it.
  
  delimiter = "\t"; //input file column delimiter
  
  ifstream infile;
    infile.open ("pressure.log");
  if(infile.fail())
  {
    cout << "Could not open pressure.log for read, terminating! \n"; //nice one!
    return 1;
  }
  while(infile.good()) // To get you all the lines, but stop on error!
  {
    getline(infile,STRING); // Saves the line in STRING.
    //sample_count++;
    lines.push_back(STRING); //did i say i love vectors?
  }
  infile.close(); //get rid of the stream
  
  lines.pop_back(); //why the heck does it put the eof-marker there?
 /* lines.pop_back();*/
  
  sample_count = lines.size(); //jep, that's the way
  //cout << "The input file vector contains " << lines.size() <<"lines!\n";
  
  
  double *samples; //let us get some variables
  fftw_complex *fftout;
  fftout = fftw_alloc_complex((sample_count)/2+1); //mem must be allocated in the first place
  samples = fftw_alloc_real(sample_count);
  
  /*for(int i = 0; i<sample_count-1; i++){
     samples[i]=0; //initialize the array with 0, for debugging purpose only 
  }*/
  
  for (int i = 0; i<lines.size(); i++){ //Hardcore C++-Verwirrung for parsing the damn input. Like perl very much, but this is horrible here...
    lines.at(i).erase(0, lines.at(i).find(delimiter) + delimiter.length());
    samples[i] = boost::lexical_cast<double> (lines.at(i));
    //cout << samples[i] << "\n";
  }
  
  cout << "Init finished \n"; //Now everything should be in place!
  
  /*for (int i = 0; i < sample_count; i++){ //Is everything there? Let's have a look at it!
    cout << i << " "<< samples[i] << "\n";
      }
      
      cout  << "Previously stands the samples array before doing anything with it!\n";*/
  
  p = fftw_plan_dft_r2c_1d(sample_count, samples, fftout,0); //There we go;
  fftw_execute(p); //Run the thing!
  fftw_destroy_plan(p); //And die!
  printFrequencies ( sample_count, (fftw_complex* )fftout ); //wanna see the results NOW!
}
