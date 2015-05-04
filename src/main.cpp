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
#include <climits>
#include <map>

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

/* Busca en la matriz la distancia entre dos puntos */
dist distancia(int i, int j, matrizDist matriz) {
    if (i < j) {
        return matriz[i][j-i-1]; 
    } else if (i == j) {
        return 0;
    } else {
        return matriz[j][i-j-1];
    }
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
    /*
        Forma de leer la matriz:
        distancia entre i y j con i < j
        matriz[i][j-i-1]
    */

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
    double pesoTotal = 0;  
    dist minDist, candidato;
  
    srand(time(NULL));

    sols.push_back(rand() % N);

    /* 
        Se guarda un mapa de distancias -> indice del punto, en el que 
        se tiene el "techo" de la probabilidad correspondiente a cada punto
        de ser el nuevo peso
    */
    
    for (int i=0; i < K-1; i++) {

        map<double, int> pesos;
        
        for (int j=0; j < N; j++) {

            minDist = LONG_MAX;
            // Obtener la distancia entre el punto j y su centro mas cercano.
            for (vector<int>::iterator it = sols.begin(); it != sols.end(); it++) {
                candidato = distancia((*it),j,matriz);
                minDist = (candidato < minDist) ? candidato : minDist;
            }

            // Guardar en el vector de pesos el nuevo punto.
            // Si la distancia es 0, quiere decir que es el mismo punto, ignorar.
            if (minDist>0) {
                // Acumular el total del peso
                pesoTotal += minDist;
                pesos[pesoTotal] = j;

            }
        }

        // Generar un numero entre 0 y 1 que representa el nuevo centro.
        // Buscar en el mapa este numero, y agregarlo a la lista de centros.
        double rnd = ((double) rand() / RAND_MAX);
        sols.push_back(pesos.upper_bound(rnd*pesoTotal)->second);

    }
    

    for (vector<int>::iterator it = sols.begin(); it != sols.end(); it++) {
        cout << (*it) << " ";
    }
    cout << endl;

    return 0;
}

