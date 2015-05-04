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
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

using namespace std;

typedef double KMcoord; //coordenada
typedef double dist;    // Distancia entre dos puntos
typedef vector< vector<dist> > matrizDist; // Matriz de distancia entre puntos.
typedef vector<KMcoord> KMpoint; // punto
typedef vector<KMpoint> KMpointArray; // arreglo de puntos

#define KM_SUM(x,y)		((x) + (y))
#define KM_POW(v)		((v)*(v))
#define K 5

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

    /*
    for(int i=0; i < dataPoints.size(); i++){
        KMpoint point = dataPoints[i];
        for(int j=0; j < point.size(); j++)
            cout << point[j] << ", " ;
        cout << endl;
                   
    }*/

    cout << "Distancia punto 98 y punto 99 -> " << kmDist(dim,
            dataPoints[98],
            dataPoints[99]) << endl;

    matrizDist matriz;

    // Llenar la matriz de distancias.
    for (int i=0; i < N-1; i++) {
        vector<dist> columna;        
        for (int j=i+1; j < N; j++) {
            columna.push_back(kmDist(dim,dataPoints[i],dataPoints[j]));
        }
        matriz.push_back(columna);
    }
    
    cout << "Distancia punto 98 y punto 99 -> " << matriz[98][0] << endl;

    // Vector de soluciones, posición en el arreglo.
    vector<int> sols;

    // Generar solucion inicial utilizando metodo k-means++
    srand(time(NULL));

    sols.push_back(rand() % N);

    for (int i=0; i < K-1; i++) {
        cout << "";
    }
    












    return 0;
}

