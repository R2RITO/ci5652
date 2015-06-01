/* 
 * File:   main.cpp
 * Authors: Artuto Voltattorni
 *          Franklin Padilla
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
#define K              3 
#define LIM_ITER       10 
#define EPSILON         0.0000005
#define DBL_MAX         1.7976931348623158e+308
#define loop(n) for(int i =0; i < n; i++) 
#define METODO_SOL_INI 0 // 0 aleatorio, 1 kmeans++

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

void imprimir_asignacion_centros(int N, pointArray dataPoints){
  cout << "Asignación de centros: " << endl;
  loop(N)
    cout << "Punto #" << i << " -> " << dataPoints[i]->centro << endl;
  
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

/*0.0.0.0
 * Dado un punto, calcula su centro más cercano.
 * @param p arregls de puntos
 *        centros vector solucion
 *        matriz Matriz de distancias
 *        N numero de puntos
 */
void calcular_centros_mas_cercanos(pointArray points, vector<int> centros,  matrizDist matriz, int N){
      
    for(int i = 0; i < N; i++ ){
        double min_d = DBL_MAX;
        for(int j = 0; j < K; j++){
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
    cout << "Método para calcular solución inicial: k-means++" << endl; 
    
    sols.push_back(rand() % N);
    /* #define loop(n) for(int ii = 0; ii < n; ++ ii)
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

void cambiar_centro(int N,
                    pointArray dataPoints,
                    vector<int> &cent,
                    matrizDist m,
                    int newc, int oldc,
                    int pos
                    )
{
  cent[pos] = newc;
  loop(N){
    if (dataPoints[i]->centro == oldc){
      dist min_d = DBL_MAX;
      for(int j = 0; j < K; j++)
      {
          dist dista = distancia(i, cent[j], m);
          if(min_d > dista)
          {
            dataPoints[i]->centro = cent[j];
            min_d = dista;
          }
      }
    } else {
      if(distancia(i, dataPoints[i]->centro, m) >= distancia(i, newc, m))
        dataPoints[i]->centro = newc;
    }
  }

}


/**
 *
 */
pair<vector<int>, dist>  calcular_mejor_vecindad(pointArray dataPoints, vector<int> sols, int N, matrizDist matriz){
    int iCenter;
    bool repetido = false, encontrado_primer_mejor = false;
    vector<int> vecino;
    double min_d = calcular_dist_solucion(dataPoints, N, matriz);
    int lim = 0;

    while(!encontrado_primer_mejor and (lim < 35)){
        vecino = sols;
        int pos_centro = rand() % (K); // elegimos un centro al azar
        while(1){
          repetido = false;
          iCenter = rand() % N;
          for (vector<int>::iterator it = sols.begin(); it != sols.end(); it++)
                if (*it == iCenter){
                  repetido = true;
                  break;
                }
          if(!repetido) break;
        }
        vecino[pos_centro] = iCenter;
        
        calcular_centros_mas_cercanos(dataPoints, vecino, matriz, N);
        double d = calcular_dist_solucion(dataPoints, N, matriz);
       
        if (d < min_d){
            min_d = d;
            encontrado_primer_mejor = true;
        }
        lim++;
    }
    pair< vector<int>, dist > result;
    result.first = (lim == 35)?sols : vecino; 
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

/**
 * Función que dado una solución, permuta sus posiciones
 *
 */
void random_shuffle(vector<int> &a){
  int r, i, temp;
  
  for(i = K-1; i > 0; i--){
    r = rand() % (i+1);
    temp = a[i];
    a[i] = a[r];
    a[r] = temp;
  }

} // fin de random_shuffle


/**
 * Local Search
 */
pair <vector <int>, dist> localSearch(int N, 
                                      pointArray dataPoints, 
                                      vector<int> sols, 
                                      matrizDist matriz)
{
// Ejecutar local search 
    int count;
    dist distorsionActual = calcular_dist_solucion(dataPoints, N, matriz);
    for (count = 0; count < LIM_ITER; count++) {
               //random_shuffle(sols); // knutt shuffle

        pair <vector<int>, dist> pr = calcular_mejor_vecindad(dataPoints, sols, N, matriz);
        vector<int> mejorVecino = pr.first;
        dist mejorVecinoDist = pr.second;

        
        if (mejorVecinoDist < distorsionActual) {
            sols = mejorVecino;
             /* Calculamos los centros más cercanos de cada punto */
            calcular_centros_mas_cercanos(dataPoints, sols, matriz, N);
            distorsionActual = mejorVecinoDist;
        }
    }
    pair< vector<int>, dist > result;
    result.first = sols; 
    result.second = distorsionActual;
    return result;

}


/**
 * Función que decide si un intercambio de centro está permitido.
 * Esto es, que el nuevo centro elegido no esté en la lista tabu y
 * que no esté repetido en la solución actual
 */
bool es_cambio_permitido(int c, vector<int> sol, vector<int> tabu){
  for(vector<int>::iterator it = sol.begin(); it != sol.end(); it++){
    if ((*it) == c)
      return false;
  }
  for(vector<int>::iterator it = tabu.begin(); it != tabu.end(); it++){
    if ((*it) == c )
      return false;
  }

  return true;
}

/**
 * Agrega un centro a la lista tabu
 */
void prohibir_centro(size_t limTabu, vector<int> &tabu, int c){
  if(tabu.size() >= limTabu)
      tabu.erase(tabu.begin());

  tabu.push_back(c);
}

/**
 * Tabú Search
 *
 */
pair< vector<int>, dist> tabuSearch(int N, 
                                    pointArray dataPoints, 
                                    vector<int> sols, 
                                    matrizDist matriz,
                                    size_t maxTabuSize){
  vector<int> tabuList, primerMejorVecino;
  vector<int> aVecino = sols;
  dist aDist;
  dist bestDist = calcular_dist_solucion(dataPoints, N, matriz);
  int c, newc, oldc, lim = 0;

  loop(LIM_ITER){
    primerMejorVecino.clear();
    for(lim = 0; lim < N/20; lim++){ // buscamos primer mejor vecino
      c = rand() % K; // elegimos un  centro al "azar"
      newc = rand() % N; // elegimos un punto al "azar"
      if (!es_cambio_permitido(newc, aVecino, tabuList)) // chequear si es tabu o repetido 
        continue;
      oldc = aVecino[c];
      cambiar_centro(N, dataPoints, aVecino, matriz, newc, oldc, c); // recalculamos centros mas cercanos
      aDist = calcular_dist_solucion(dataPoints, N, matriz); // distorsion de la nueva solucion

      if (aDist < bestDist)
      {
        primerMejorVecino = aVecino;
        bestDist = aDist;
        prohibir_centro(maxTabuSize, tabuList, oldc);
        break; // encontrado mejor Vecino
      } 
      else 
      {
        prohibir_centro(maxTabuSize, tabuList, newc);
        cambiar_centro(N, dataPoints, aVecino, matriz, oldc, newc, c); 
      }

      aVecino = sols;

    }

    if (lim != N/20){
      sols = primerMejorVecino;
    } else {
      /* calcular_centros_mas_cercanos(dataPoints, sols, matriz, N); */
    }

    lim = 0;
    
  }
  pair< vector<int>, dist> result;
  result.first = sols;
  result.second = bestDist;
  return result;
}


dist dbIndex(int N, pointArray dataPoints, vector<int> centers, matrizDist m){
  vector<dist> S; // distorsion de cada cluster
  vector<dist> D;
  vector<int> counters; // para guardar el tamaño de cada cluster
  map<int, int> centIndex; // mapea posicion de centros dentro del vector de solucion
  loop(K){ 
    centIndex[centers[i]] = i; 
    S.push_back(0.0); 
    counters.push_back(0);
    D.push_back(0.0);
  }

  loop(N){ // para cada centro, sumamos la distancias de sus puntos mas cercanos
    int centro = dataPoints[i]->centro;
    int index = centIndex[centro];
    S[index] = S[index] + distancia(i, centro, m); 
    counters[index] = counters[index] + 1; // aumentamos los contadores de tamaño
  }
  
  loop(K)
    S[i] = S[i] / counters[i]; // dividimos entre tamaño de cluster

  for(int i = 0; i < K-1; i++){
    for(int j = i+1; j < K;  j++){
      D[i] = max(D[i], (S[i] + S[j]) / distancia(centers[i], centers[j], m));
    }
  }

  dist DB = 0.0;
  loop(N)
    DB += D[i];

  return DB/N;
  
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
    srand(time(NULL));

    //cálculo solución inicial
    if (METODO_SOL_INI)
       sols = kpp(N, matriz); // kmeans++ 
    else {
      /* cout << "Método para calcular solución inicial: aleatorio" << endl; */ 
      loop(K)
        sols.push_back(rand() % N);
    }

    calcular_centros_mas_cercanos(dataPoints, sols, matriz, N); // calculamos centros mas cercanos
    dist distInicial = calcular_dist_solucion(dataPoints, N, matriz);
    /* cout << "Nro Iteraciones en el LS: " << LIM_ITER << endl; */
    /* cout << "Solucion inicial: " ; */
    /* loop(K) */
    /*   cout << sols[i] << "; "; */
    /* cout << endl << "Distorsion inicial: " << distInicial << endl; */
    /* cout << "Davis Goulding Index de la solucion inicial: " << dbIndex(N, dataPoints, sols, matriz) << endl; */
    dist dbiInicial = dbIndex(N, dataPoints, sols, matriz); 

    /* pair <vector<int>, dist> pr = localSearch(N, dataPoints, sols, matriz); //aplicamos LocalSearch */
    pair <vector<int>, dist> pr = tabuSearch(N, dataPoints, sols, matriz, N/4); // o aplicamos tabuSearch 
    sols = pr.first;
    dist distFinal = pr.second;

    /* cout << "Solucion Final: "; */
    /* loop(K) */
    /*   cout << sols[i] <<  "; " ; */
    /* cout << endl << "Distorsión Final: " << distFinal << endl; */
    /* cout << "Davis Goulding Index de la solucion final: " << dbIndex(N, dataPoints, sols, matriz) << endl; */
    /* cout << "Porcentaje de mejora con respecto a la sol inicial: " << (distInicial-distFinal)/distInicial << endl; */
   
    calcular_centros_mas_cercanos(dataPoints, sols, matriz, N); 
    
    dist dbiFinal =  dbIndex(N, dataPoints, sols, matriz);
    cout << distFinal << " " << (distInicial-distFinal)/distInicial << " ";
    cout << dbiFinal << " " << (dbiInicial-dbiFinal) / dbiInicial << endl; 

    return 0;
}
