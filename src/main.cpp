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
#define K               5
#define LIM_ITER        25
#define EPSILON         0.0000005
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


 
pair<vector<int>, dist>  calcular_mejor_vecindad(pointArray dataPoints, vector<int> sols, int N, matrizDist matriz){
    int iCenter;
    bool repetido;
    vector<int> mejorVecino = sols;
    srand(time(NULL));
    double min_d = calcular_dist_solucion(dataPoints, N, matriz);
    
    for(int i = 0; i < K; i++){
        vector<int> vecino;
        vecino = sols;
        repetido = true;
        while(repetido){
            iCenter = rand() % N;
            for (vector<int>::iterator it = sols.begin(); it != sols.end(); it++)
                if (*it == iCenter)
                    continue;
            repetido = false;
        }
        vecino[i] = iCenter;
        
        calcular_centros_mas_cercanos(dataPoints, vecino, matriz, N);
        double d = calcular_dist_solucion(dataPoints, N, matriz);
        if (d < min_d){
            min_d = d;
            mejorVecino = vecino;           
        }
    }
    pair< vector<int>, dist > result;
    result.first = mejorVecino;
    result.second = min_d;
    return result;
}

/* Determina si hubo convergencia */
int convergencia(dist d1, dist d2) {

    if (labs(d1-d2) < EPSILON) {
        return 1;
    } else {
        return 0;
    }

}


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
        stringstream lineStream(lineData);
        while(getline(lineStream, subLineData, ',')){ // Dividimos por comas
            coord num = atof(subLineData.c_str());
            p->coords.push_back(num); 
        }
        dataPoints.push_back(p);
    }
    
    N = dataPoints.size();
    dim = (dataPoints[0]->coords).size();

    matrizDist matriz;

    // Llenar la matriz de distancias.
    for (int i=0; i < N-1; i++) {
        vector<dist> columna;        
        for (int j=i+1; j < N; j++) {
            columna.push_back(kmDist(dim,dataPoints[i],dataPoints[j]));
        }
        matriz.push_back(columna);
    }

    // Vector de soluciones, posición de cada centro en el arreglo de puntos.
    vector<int> sols;
    sols.clear();
    
    // Generar solucion inicial utilizando metodo k-means++
    sols = kpp(N, matriz); 
    calcular_centros_mas_cercanos(dataPoints, sols, matriz, N);
    dist distorsionActual = calcular_dist_solucion(dataPoints, N, matriz);
    dist mejorVecinoDist = 0;

    // Ejecutar local search 
    int count;
    for (count = 0; count < LIM_ITER; count++) {
        /* Calculamos los centros más cercanos de cada punto */
        calcular_centros_mas_cercanos(dataPoints, sols, matriz, N);         
        pair <vector<int>, dist> pr = calcular_mejor_vecindad(dataPoints, sols, N, matriz);
        vector<int> mejorVecino = pr.first;
        dist mejorVecinoDist = pr.second;

        if (convergencia(distorsionActual,mejorVecinoDist))
            break;

        if (mejorVecinoDist < distorsionActual) {
            sols = mejorVecino;
            distorsionActual = mejorVecinoDist;
        }
    }

    cout << "Solucion: ";
    for (int i=0; i < K; i++) 
        cout << sols[i] <<  "; " ;
    cout << endl;
    
    return 0;
}
