#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <map>
#include <sys/time.h>

#define NUM_ITER 15

using namespace std;

int main(int argc, char* argv[]) {

    struct timeval diff, startTV, endTV;
    long totalSec = 0;
    long totalMSec = 0;

    if(argc < 2){
      cerr << "Argumento invalido" << endl;
      return 1;
    }

    string exec = "./main < ";
    string commandStr = exec + argv[1];

    for (int i = 0; i < NUM_ITER; i++) {

        gettimeofday(&startTV, NULL); 

        std::system(commandStr.c_str());

        gettimeofday(&endTV, NULL); 

        timersub(&endTV, &startTV, &diff);

        totalSec += diff.tv_sec;
        //cout << diff.tv_usec << endl;
        totalMSec += diff.tv_usec;
    }
    
    //cout << totalMSec << endl;
    //printf("time taken = %ld %ld\n", totalSec, totalMSec);
    
    //cout << endl;

}
