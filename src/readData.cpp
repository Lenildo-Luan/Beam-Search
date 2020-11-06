#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;



void readData(int argc, char** argv, int *n, double ***mJobs, double ***mSetupTimes){

    if (argc < 3) {
         cout << "\nFaltando parametros\n";
         cout << " ./exec [Instancia] "<< endl;
         exit(1);
     }

     if (argc > 4) { // testes insertion
          cout << "\nMuitos parametros\n";
          cout << " ./exec [Instancia] " << endl;
         exit(1);
     }

    int N, nMachines;
    string arquivo, ewt;

    char *instancia;
    instancia = argv[1];
    int lixo;


    ifstream in( instancia, ios::in);

    in >> N;
    in >> nMachines;
    for(int i = 0; i < 2; i++){
        in >> lixo;
    }
    *n = N;
    int mteste[N][4];
    //int initialSetupArray[N];

    double **jobsInf = new double*[N+1];

    for ( int i = 1; i <= N; i++ ) {
        jobsInf[i] = new double [4];
    }

   // int *initialSetupArray = new int[N+1];

    double **setupTimes = new double*[N+1];

    for ( int i = 0; i <= N; i++){
        setupTimes [i] = new double [N+1];
    }

    //int setupTimes[N][N];


    for(int i = 1; i <= N; i++){
        for(int j = 0; j < 4; j++){
            in >> jobsInf[i][j];
        }
    }

    for(int i = 0; i <= N; i++){
        setupTimes[i][0] = 0; // sem erro, ainda suspeito
        for(int j = 1; j <= N; j++){
            in >> setupTimes[i][j];
        }
    }

    /*for(int i = 0; i <= N; i++){
        setupTimes[i][0] = 0;
        for(int j = 1; j <= N; j++){
            setupTimes[i][j] = 0;
        }
    }*/



   *mSetupTimes = setupTimes;
   *mJobs = jobsInf;

}