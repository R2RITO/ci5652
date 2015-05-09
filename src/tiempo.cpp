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

#define NUM_ITER 1 

using namespace std;

int main() {

    struct timeval diff, startTV, endTV;
    long totalSec = 0;
    long totalMSec = 0;

    for (int i = 0; i < NUM_ITER; i++) {

        gettimeofday(&startTV, NULL); 

        std::system("./main < iris.data.txt");

        gettimeofday(&endTV, NULL); 

        timersub(&endTV, &startTV, &diff);

        totalSec += diff.tv_sec;
        totalMSec += diff.tv_usec;
    }
    
    cout << endl;
    printf("time taken = %ld %ld\n", totalSec, totalMSec);
    
    cout << endl;

}
