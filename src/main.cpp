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
#include <map>

using namespace std;

typedef double coord; //coordenada
typedef double dist;    // Distancia entre dos puntos
typedef vector< vector<dist> > matrizDist; // Matriz de distancia entre puntos.
typedef struct { vector<coord> coords; int centro; } point_t, *point; // punto
typedef vector<point> pointArray; // arreglo de puntos

#define KM_SUM(x,y)		((x) + (y))
#define KM_POW(v)		((v)*(v))
#define K 5
#define DBL_MAX         1.7976931348623158e+308

/* 
 * Calcula la distancia euclideana de 2 puntos
 * 
 */
dist kmDist(			// interpoint squared distance
    int			dim,
    point		p,
    point		q)
{
    int d;
    double diff;
    dist distance = 0;

    for (d = 0; d < dim; d++) {
	diff = (p->coords)[d] - (q->coords)[d];
	distance = KM_SUM(distance, KM_POW(diff));
    }
    return distance;
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
 * Dado un punto, calcula su centro más cercano.
 * @param p arregls de puntos
 *        centros vector solucion
 *        matriz Matriz de distancias
 *        N numero de puntos
 */
void calcular_centros_mas_cercanos(pointArray points, vector<int> centros,  matrizDist matriz, int N){
      
    for(int i = 0; i < N; i++ ){
        double min_d = DBL_MAX;
        for(int j = 0; j < K-1; j++){
            dist dista = distancia(i, centros[j], matriz);
            if(min_d > dista){
                points[i]->centro = centros[j];
                min_d = dista;
            }
        }
    }
}

/**
 * Funcion que genera una solucion inicial usando el metodo de k-means
 * @param N numero de puntos
 * @param sols vector solucion
 * @param matriz 
 */
vector<int> kpp(int N, matrizDist matriz){
 
    double pesoTotal = 0;  
    dist minDist, candidato;
    vector<int> sols;
    
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

            minDist = DBL_MAX;
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
    return sols;
}

/**
 * Calcula la distorsion de una solucion
 * @param dataPoints
 * @param N
 * @param matriz
 * @return 
 */
dist calcular_dist_solucion(pointArray dataPoints, int N, matrizDist matriz){
    dist res = 0;
    for(int i = 0; i < N; i++){
        res += distancia(i, dataPoints[i]->centro, matriz);
    }
    return res;
}

/**
 * 
 
vector< vector<int> > calcular_vecindad(vector<int> sols){
    vector< vector<int> >vecindad;
    for(int i = 0; i < K; i++){
        vector<int> vecino;
        vector
    }

}
*/

/*
 * 
 */
int main(int argc, char** argv) {
    
   
    int dim = 0, N = 0;
    string lineData; 
    string subLineData;
    pointArray dataPoints; 
    
    while( getline(cin, lineData)) {
        point p = new point_t();
        p->coords.clear();
        cout << "size:" << p->coords.size();
        stringstream lineStream(lineData);
        while(getline(lineStream, subLineData, ',')){ // Dividimos por comas
            coord num = atof(subLineData.c_str());
            p->coords.push_back(num); 
        }
        dataPoints.push_back(p);
    }
    
    N = dataPoints.size();
    dim = (dataPoints[0]->coords).size();
    
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

    // Vector de soluciones, posición de cada centro en el arreglo de puntos.
    vector<int> sols;
    sols.clear();

    /* Elegimos aleatoriamente la solucion inicial */
    //srand(time(NULL)); 
    //for (int i=0; i < K-1; i++) 
     //   sols.push_back(rand() % N);
    
    
    // Generar solucion inicial utilizando metodo k-means++
    sols = kpp(N, matriz);   
    
    /* Calculamos los centros más cercanos de cada punto */
    calcular_centros_mas_cercanos(dataPoints, sols, matriz, N);
    
    cout << "solucion inicial: ";
    for (int i=0; i < K-1; i++) 
        cout << sols[i] <<  "; " ;
    
    cout << endl;
    for (int i=0; i < N; i++) 
        cout << "punto n°" << i << " centro: " << dataPoints[i]->centro << endl;
    
    cout << "distorsion solucion: " << calcular_dist_solucion(dataPoints, N, matriz) << endl;
    
    return 0;
}
