/* 
 * File:   main.cpp
 * Author: franklin
 *
 * Created on May 2, 2015, 2:44 PM
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

typedef double KMcoord; //coordenada
typedef vector<KMcoord> KMpoint; // punto
typedef vector<KMpoint> KMpointArray; // arreglo de puntos

#define KM_SUM(x,y)		((x) + (y))
#define KM_POW(v)		((v)*(v))

/* 
 * Calcula la distancia euclideana de 2 puntos
 * 
 */
double kmDist(			// interpoint squared distance
    int			dim,
    KMpoint		p,
    KMpoint		q)
{
    int d;
    KMcoord diff;
    KMcoord dist;

    dist = 0;
    for (d = 0; d < dim; d++) {
	diff = p[d] - q[d];
	dist = KM_SUM(dist, KM_POW(diff));
    }
    return dist;
}

/*
 * 
 */
int main(int argc, char** argv) {
    
   
    int dim = 0, N = 0;
    string lineData; 
    string subLineData;
    KMpointArray dataPoints; 
    
    while( getline(cin, lineData)) {
        KMpoint p;
        stringstream lineStream(lineData);
        while(getline(lineStream, subLineData, ',')){ // Dividimos por commas
            KMcoord num = atof(subLineData.c_str());
            p.push_back(num);          
        }
        dataPoints.push_back(p);
    }
    
    N = dataPoints.size();
    dim = dataPoints[0].size();
    
    cout << "Numero de datos: " << N << endl << "Dimensiones: " << dim << endl;
    
    for(int i=0; i < dataPoints.size(); i++){
        KMpoint point = dataPoints[i];
        for(int j=0; j < point.size(); j++)
            cout << point[j] << ", " ;
        cout << endl;
                   
    }
    cout << "Distancia punto 1 y punto 2 -> " << kmDist(dim,
            dataPoints[1],
            dataPoints[2]);
    
    
    return 0;
}

