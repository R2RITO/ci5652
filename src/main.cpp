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
#include <math.h> 

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
#define DBL_MIN        -DBL_MAX
#define loop(n) for(int i =0; i < n; i++) 
#define METODO_SOL_INI 0 // 0 aleatorio, 1 kmeans++
#define TAM_POBLACION       16
#define MAX_GENERACIONES    15
#define PROB_CRUCE          0.4
#define PROB_MUTACION       0.12
#define TAM_TORNEO          3

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
    //cout << "Método para calcular solución inicial: k-means++" << endl; 
    
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


void sumarPuntos(vector< vector<coord> > *p, point q, int indx) {
    int dim = q->coords.size();

    for (int i=0; i < dim; i++) {
        (*p)[indx][i] = (*p)[indx][i] + q->coords[i];
    }

    return;
}

/**
 * Calcula las medias de los clusters (lloyd)
 * @param dataPoints
 * @param K
 * @param sols
 * @return vector con las medias.
 */
vector< vector<coord> > calcularMedias(pointArray dataPoints, vector< vector<coord> > sols, int N) {

    // Inicializar las medias.
    vector< vector<coord> > vectorMedias;
    vector<int> numElems;

    for (int i=0; i < K; i++) {
        vector<coord> comp;
        coord dim = 0.0;
        for (size_t j=0; j < sols[0].size(); j++) {
            comp.push_back(dim);
        }
        vectorMedias.push_back(comp);
        numElems.push_back(0);
    }

    for (int i=0; i < N; i++) {
        sumarPuntos(&vectorMedias, dataPoints[i], dataPoints[i]->centro);
        numElems[dataPoints[i]->centro] += 1;
    }

    for (int i=0; i < K; i++) {
        for (size_t j=0; j<vectorMedias[0].size(); j++) {
            vectorMedias[i][j] = vectorMedias[i][j] / numElems[i];
        }
    }

    return vectorMedias;

}

/**
 * Calcula el centro correspondiente para cada punto, retorna la distorsión
 * total.
*/

dist calcularCentrosCercanosLloyd(pointArray dataPoints, vector< vector<coord> > sols, int N) {

    // Crear point para poder utilizar la funcion kmDist.
    pointArray solsP;    
    for (int h = 0; h < K; h++) {
        point p = new point_t();
        p->coords = sols[h];
        p->centro = 0;
        solsP.push_back(p);
    }

    int dim = sols[0].size();
    dist distorsionTotal = 0.0;

    //cout << "Centro de 0 antes: " << dataPoints[0]->centro << endl;

    for(int i = 0; i < N; i++ ){
        double min_d = DBL_MAX;
        for(int j = 0; j < K; j++){
            dist dista = kmDist(dim, dataPoints[i], solsP[j]);
            if(min_d > dista){
                dataPoints[i]->centro = j;
                min_d = dista;
            }
        }
        distorsionTotal += min_d;
    }

    //cout << "Centro de 0 después: " << dataPoints[0]->centro << endl;

    return distorsionTotal;
}

/**
 * Calcula una solucion utilizando el algoritmo de Lloyd
 * @param dataPoints
 * @param N
 * @param matriz
 * @return vector con la solucion, distorsion de la solucion.
 */
pair<vector< vector<coord> >, double> lloyd(vector< vector<coord> > solIni, pointArray dataPoints, int N) {

    dist distActual;
    vector< vector<coord> > sols = solIni;

    for (int count = 0; count < LIM_ITER; count++) {
        // Calcular nuevos centroides.
        sols = calcularMedias(dataPoints, sols, N);

        // Reasignar cada punto.
        distActual = calcularCentrosCercanosLloyd(dataPoints, sols, N);
        //cout << "Dist actual de lloyd: " << distActual << endl;        

    }

    pair< vector< vector<coord> >, dist > result;
    result.first = sols; 
    result.second = distActual;
    return result;
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


/* Funcion que utiliza kpp para calcular una solucion inicial y la adapta al formato
 * utilizado por lloyd para calcular la solucion final.
*/
pair<vector< vector<coord> >, double> solucionLloyd(pointArray dataPoints, int N, matrizDist matriz, vector<int> sols) {

    /* Adaptar los centros, y asignar a cada punto su centro */
    vector< vector<coord> > solsAdapt;
    for (int i = 0; i < K; i++) {
        vector<coord> centro = dataPoints[sols[i]]->coords;
        solsAdapt.push_back(centro);
    }

    dist distIni = calcularCentrosCercanosLloyd(dataPoints, solsAdapt, N);

    cout << "Dist init de lloyd " << distIni << endl;

    return lloyd(solsAdapt, dataPoints, N);

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

/* 
################################################################
################################################################
################################################################
            Seccion de metaheuristica GRASP
################################################################
################################################################
################################################################
*/

/* Construye una solucion candidata utilizando la estrategia de grasp */

vector<int> construirSolGrasp(double alpha, matrizDist matriz, int N) {

    // Seleccionar el primer elemento aleatorio
    vector<int> sols;    
    sols.push_back(rand() % N);

    pair<int, double> cMin, cMax, actual;
    vector< pair<int, double> > RCL;
    dist minDist, maxRCL, minRCL, distCandidato;

    for (int i = 1; i < K; i++) {

        /* Reiniciar RCL. maximo y minimo */
        RCL.clear();
        cMax.first = -1;
        cMax.second = DBL_MIN;
        cMin.first = -1;
        cMin.second = DBL_MAX;
    
        /*  Calcular la distancia de cada elemento a los centroides, guardarlo
            en el RCL.
        */
        
        for (int j=0; j < N; j++) {

            minDist = DBL_MAX;
            // Obtener la distancia entre el punto j y su centro mas cercano.
            for (vector<int>::iterator it = sols.begin(); it != sols.end(); it++) {
                distCandidato = distancia((*it),j,matriz);
                minDist = (distCandidato < minDist) ? distCandidato : minDist;
            }

            // Si la distancia es 0, quiere decir que es el mismo punto, ignorar.
            if (minDist > 0) {
                actual.first = j;
                actual.second = minDist;

                /* Actualizar el minimo y el maximo */
                if (minDist < cMin.second) {
                    cMin.first = j;
                    cMin.second = minDist;
                }

                if (minDist > cMax.second) {
                    cMax.first = j;
                    cMax.second = minDist;
                }

                /* Guardar el elemento en el RCL */
                RCL.push_back(actual);
            }
        }

        /* Seleccionar un elemento aleatorio del RCL */
        maxRCL = cMax.second;
        minRCL = cMax.second - alpha*(cMax.second - cMin.second);
        int centro;
        bool encontrado = false;

        /* Seleccionar un elemento aleatorio del RCL, si cumple con
           la condicion de bondad (su distancia al centroide mas cercano
           es mayor a minRCL, entonces sera el centro a agregar.
        */

        while (!encontrado) {
            centro = rand() % (N - i);
            encontrado = (RCL[centro].second > minRCL) ? true : false;
        }

        sols.push_back(centro);

    }

    return sols;

}


pair<vector<int>, dist> grasp(matrizDist matriz, int N, pointArray dataPoints, int lim) {

    pair<vector<int>, dist> solCandidata, optimo;
    optimo.second = DBL_MAX;

    double alpha;    
    
    // Crear el vector de probabilidades para alpha
    map<double, double> valoresAlpha;
    double acc = 0.1;
    for (int j = 0; j < 9; j++) {
        valoresAlpha[acc] = acc;
        acc += 0.1;
    }
     

    for (int i = 0; i < lim; i++) {

        // Seleccion de alpha.
        alpha = valoresAlpha.upper_bound(rand() / RAND_MAX)->second;

        // Fase de construccion
        vector<int> candidato = construirSolGrasp(alpha, matriz, N);

        // Fase de mejora
        calcular_centros_mas_cercanos(dataPoints, candidato, matriz, N);
        solCandidata = localSearch(N,dataPoints,candidato,matriz);

        // Actualizar optimo
        if (solCandidata.second < optimo.second) {
            optimo.first = solCandidata.first;
            optimo.second = solCandidata.second;
        }  

        //cout << "Voy por la i: " << i << endl; 

    }

    return optimo;

}

/* Fin de GRASP */

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

/******************************************
**                                       **
**              Genetico                 **
**                                       **
******************************************/ 

double calcularDistorsionGenetico(pointArray points, vector<int> centros,  matrizDist matriz, int N) {

    double disTotal = 0.0;

    for(int i = 0; i < N; i++ ){
        dist min_d = DBL_MAX;
        dist dista;
        for(int j = 0; j < K; j++){
            dista = distancia(i, centros[j], matriz);
            min_d = (dista < min_d) ? dista : min_d;
        }

        disTotal += min_d;
    }

    return disTotal;
}

int seleccionTorneo(vector< vector<int> > pob, matrizDist matriz, pointArray dataPoints, int T, int N) {

    // Generar los "contendientes" del torneo
    // El par contiene el indice de la solucion de la poblacion y su distorsion
    pair<int, double> distorsionMax, distorsionActual;
    distorsionMax.second = DBL_MAX;
    int contendiente;

    for (int i = 0; i < T; i++) {

        contendiente = rand() % TAM_POBLACION;
        distorsionActual.first = contendiente;
        distorsionActual.second = calcularDistorsionGenetico(dataPoints, pob[contendiente], matriz, N);
        distorsionMax = (distorsionActual.second < distorsionMax.second) ? distorsionActual : distorsionMax;

    }

    return distorsionMax.first;
}

pair< vector<int>, double > genetico_generacional(matrizDist matriz, int N, pointArray dataPoints) {
	
    // Mapa de probabilidades de cruce, fijo para todas las iteraciones.
    map< double, int > mapCruces;
    vector<int> sol1 (K,0);
    vector<int> sol2 (K,0);
    vector< vector<int> > nuevaPob;
    double acc = 0.0;
    int padA, padB;
    double rnd, aux, minDist, distAct;
    vector<int> mascara (K,0);
    vector<int> mutacion (K,0);
    vector<int> solIni, cruzandos;

    vector< vector<int> > poblacion;

    for ( int s = 0; s < TAM_POBLACION; s++) {
        acc += PROB_CRUCE;
        mapCruces[acc] = s;

        // Seccion para generar poblacion inicial
        solIni.clear();
        for (int sIndx = 0; sIndx < K; sIndx++) {
            solIni.push_back(rand() % N);
        }

        poblacion.push_back(solIni);

    }

    // Inicio del algoritmo

    for( int numGens = 0; numGens < MAX_GENERACIONES; numGens++ ) {

        // Seleccionar
        // Se utiliza la estrategia de torneo.
        // O(n*k*tamTorneo) para cada padre. Esto ejecuta pob veces.

        // cruzandos contiene los indices de las soluciones en la poblacion.
        cruzandos.clear();

        for ( int pad = 0; pad < TAM_POBLACION; pad++ ) {
            cruzandos.push_back(seleccionTorneo(poblacion, matriz, dataPoints, TAM_TORNEO, N));
        }

        // Cruzar

        // Generar la máscara que servirá para el cruce uniforme.
        for ( int i = 0; i < K; i++ ) {
            rnd = ((double) rand() / RAND_MAX);
            mascara[i] = (rnd < 0.5) ? 0xFFFF : 0x0000;
            mutacion[i] = (rnd < 0.35) ? 0xFFFF : 0x0000;
        }

        nuevaPob.clear();

        for ( int i = 0; i < TAM_POBLACION / 2; i++ ) {

            // Selecciono dos padres aleatorios para cruzar
            rnd = ((double) rand() / RAND_MAX);
            padA = mapCruces.upper_bound(rnd*acc)->second;

            rnd = ((double) rand() / RAND_MAX);
            padB = mapCruces.upper_bound(rnd*acc)->second;

            // Cruzo y genero dos nuevas soluciones, que reemplazan a sus padres.
            // Iterar para aplicar las máscaras a cada componente de las soluciones.

            for (int h = 0; h < K; h++ ) {
                sol1[h] = (poblacion[cruzandos[i]][h] & mascara[h]) + (poblacion[cruzandos[i + (TAM_POBLACION / 2)]][h] & ~mascara[h]);
                sol2[h] = (poblacion[cruzandos[i]][h] & ~mascara[h]) + (poblacion[cruzandos[i + (TAM_POBLACION / 2)]][h] & mascara[h]);
            }

            nuevaPob.push_back(sol1);
            nuevaPob.push_back(sol2);
            
        }

        // Mutar
        // Elijo el numero de sols a mutar, para cada una selecciono una sol
        // aleatoria y con prob 0.35 muto sus componentes a otra aleatoria.
        for ( int nmut = 0; nmut < TAM_POBLACION * PROB_MUTACION; nmut++ ) {
            rnd = (rand() % TAM_POBLACION);
            
            for ( int mut = 0; mut < K; mut++ ) {
                aux = ((double) rand() / RAND_MAX);
                sol1[mut] = (aux < 0.35) ? floor(aux * (N-1)) : nuevaPob[rnd][mut];                         
            }

            nuevaPob[rnd] = sol1;
        }

        // Aumentar probabilidad de mutacion
        // Disminuir probabilidad de cruce

        poblacion = nuevaPob;
    }

    pair< vector<int>, double > mejorSol;
    mejorSol.second = DBL_MAX;

    for (int m = 0; m < TAM_POBLACION; m++ ) {

        distAct = calcularDistorsionGenetico(dataPoints, poblacion[m], matriz, N);
        
        if (distAct < mejorSol.second) {
            mejorSol.first = poblacion[m];
            mejorSol.second = distAct;
        } 
    }

    return mejorSol;


}

/* Fin de genetico */

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

    vector<int> solsInit = sols;

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
    /*pair <vector<int>, dist> pr = tabuSearch(N, dataPoints, sols, matriz, N/4); // o aplicamos tabuSearch 
    sols = pr.first;
    dist distFinal = pr.second;*/

    /* cout << "Solucion Final: "; */
    /* loop(K) */
    /*   cout << sols[i] <<  "; " ; */
    /* cout << endl << "Distorsión Final: " << distFinal << endl; */
    /* cout << "Davis Goulding Index de la solucion final: " << dbIndex(N, dataPoints, sols, matriz) << endl; */
    /* cout << "Porcentaje de mejora con respecto a la sol inicial: " << (distInicial-distFinal)/distInicial << endl; */
   
    /*
    calcular_centros_mas_cercanos(dataPoints, sols, matriz, N); 
    
    dist dbiFinal =  dbIndex(N, dataPoints, sols, matriz);
    cout << distFinal << " " << (distInicial-distFinal)/distInicial << " ";
    cout << dbiFinal << " " << (dbiInicial-dbiFinal) / dbiInicial << endl; 
    */

    /*
    pair< vector< vector<coord> >, dist > resultLloyd = solucionLloyd(dataPoints, N, matriz, solsInit);
    dist lloydFinal = resultLloyd.second;
    cout << endl << "Distorsión Final de Lloyd: " << lloydFinal << endl;
    */

    /*
    pair <vector<int>, dist> gp = grasp(matriz, N, dataPoints, 6);
    //cout << "Distorsión final de GRASP: " << gp.second << endl;
    //cout << "Porcentaje de mejora con respecto a la sol inicial: " << (distInicial-gp.second)/distInicial << endl;
    cout << gp.second << " " << (distInicial-gp.second)/distInicial << " ";
    */

    pair <vector<int>, dist> gen = genetico_generacional(matriz, N, dataPoints);
    //cout << "Distorsión final de Genético: " << gen.second << endl;
    cout << gen.second << " ";
    //cout << "Porcentaje de mejora con respecto a la sol inicial: " << (distInicial-gp.second)/distInicial << endl;
    //cout << gp.second << " " << (distInicial-gp.second)/distInicial << " ";

    return 0;
}
