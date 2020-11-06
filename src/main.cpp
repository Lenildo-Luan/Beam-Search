#include "readData.h"
#include "CustoIn.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include "infoSeq.h"
#include <sys/timeb.h>
#include <sys/resource.h>
#include <fstream>
#include <tuple>

using namespace std;

ofstream out("arquivo.txt");
double ** mJobs; // matriz de adjacencia
double ** mSetupTimes; // matriz reorganizada;
int n; // quantidade total de vertices
vector<double> compTimes;

int melhoras = 0, melhorasSwap = 0, melhoras2opt = 0;
int melhorasptb = 0;
int melhorasReinsert[13];
vector<int> listaSub = {0,1,2,3};

//mJobs[1][] = release date
//mJobs[2][] = processing time

void printData(int n, double ** mJobs, double ** mSetupTimes);
vector<int> construction(int n, double ** mJobs, double ** mSetupTimes, double alfa, double &cost);

//Esse aqui
double sequenceTime(vector<int> s, vector<int> remainingVertices, double ** mJobs, double **mSetupTimes);
double idleTime(vector<int> s, double cTime, int node, double ** mJobs, double **mSetupTimes);
//Esse aqui

vector<CustoIn> calculaCusto(vector<int> listaCandidatos, vector<int> &s, double custoAtual);
bool comp(const CustoIn& a, const CustoIn& b);
double compCostSwap(vector<int> &s, int i, int j, double &compTime);
//void swap(vector<int> &solution, double &cost);
double compCostBlock(vector<int> &s, int lb, int a, int b, double compTime);
void compCompletionTime(vector<int> &s, vector<double> &compTimes);
vector<int> perturb(vector<int> solucao);
vector<int> ils(int i_max, int i_ils);
void reInsertion(int l, vector<int> &solucao, double &custo);
double compCostReinsertion(int l, int n, vector<int> &s, int i, int j, vector<double> compTimes);
void rvnd(vector<int> &solucao, double &custo);
double revCompCostBlock(vector<int> &s, int lb, int a, int b, double compTime);
void twoOptN(vector<int> &solucao, double &custo);
double compCost2opt(vector<int> &s, int i, int j, vector<double> compTimes);
double compCostReinsertionv2(int l, vector<int> &s, int i, int j, infoSequence **sequencesMatrix);
double cpuTime();
void swap(vector<int> &solution, double &cost, infoSequence **sequencesMatrix);
void printSolution(vector<int> solucao, double **mJobs, double ** mSetupTimes);

//Luan implementations
void beamSearchLevel(pair< vector<int>, double> &bestSol, int totalBranchesPerNode, int totalBranchesPerLevel);
void beamSearchDepth(pair< vector<int>, double> &bestSol, pair< vector<int>, double> sol, vector<int> remainigVertices, 
                vector<int> &branchesPerLevel, int totalBranchesPerNode, int totalBranchesPerLevel, int currentLevel);
int binarySearch(vector<int> &vec, int posInit, int posEnd, int compare);
bool compBS(pair< vector<int>, double > a, pair< vector<int>, double > b);
bool compBSLWT(tuple <pair<double, double>, vector<int>, vector<int> > a, tuple <pair<double, double>, vector<int>, vector<int> > b);
bool compBSLCT(tuple <pair<double, double>, vector<int>, vector<int> > a, tuple <pair<double, double>, vector<int>, vector<int> > b);
//End

int main(int argc, char** argv) {

    readData(argc, argv, &n, &mJobs, &mSetupTimes);

    vector<int> remainigVertices, branchesPerLevel;
    pair< vector<int>, double> bestSol, sol;
    bestSol.second = numeric_limits<double>::max();
    sol.first = {};
    sol.second = 0;

    for (int i = 1; i <= n; i++){
        remainigVertices.push_back(i);
        branchesPerLevel.push_back(0);
    }

    // vector<int> test = {2, 6, 4, 3}, test2 = {1, 5};
    // sequenceTime(test, test2, mJobs, mSetupTimes);

    beamSearchLevel(bestSol, 5, 500);
    printSolution(bestSol.first, mJobs, mSetupTimes);

    out.close();

    return 0;
}

//Luan implementations
void beamSearchLevel(pair< vector<int>, double> &bestSol, int totalBranchesPerNode, int totalBranchesPerLevel){
    int removePosition = 0;

    vector<int> listOfNodesRemaining;
    for (int i = 1; i <= n; i++) listOfNodesRemaining.push_back(i);

    vector <tuple<pair<double, double>, vector<int>, vector<int> >> listOfSols;
    vector <tuple<pair<double, double>, vector<int>, vector<int> >> nodeListOfSols;
    vector <tuple<pair<double, double>, vector<int>, vector<int> >> levelListOfSols;

    vector <tuple<pair<double, double>, vector<int>, vector<int> >> emptSols;
    tuple <pair<double, double>, vector<int>, vector<int> > sol({0, 0}, {}, listOfNodesRemaining);
    listOfSols.push_back(sol);                                                  

    for (int i = 0; i < n; i++){
        cout << "Nivel: " << i << endl;
        levelListOfSols.clear();
        for (int j = 0; j < listOfSols.size(); j++){
            // cout << nodeListOfSols.size() << endl;
            nodeListOfSols.clear();
            for (int k = 0; k < get<2>(listOfSols[j]).size(); k++){
                nodeListOfSols.push_back(sol); 
                nodeListOfSols[k] = listOfSols[j];
                // get<2>(nodeListOfSols[k]) = get<2>(listOfSols[j]);
                // get<1>(nodeListOfSols[k]) = get<1>(listOfSols[j]);
                get<0>(nodeListOfSols[k]).first = idleTime(get<1>(nodeListOfSols[k]), get<0>(nodeListOfSols[k]).second, 
                                                           get<2>(listOfSols[j])[k], mJobs, mSetupTimes);
                get<1>(nodeListOfSols[k]).push_back(get<2>(listOfSols[j])[k]);
                removePosition = binarySearch(get<2>(nodeListOfSols[k]), 0, get<2>(nodeListOfSols[k]).size()-1,
                                              get<1>(nodeListOfSols[k])[ get<1>(nodeListOfSols[k]).size()-1 ]);
                get<2>(nodeListOfSols[k]).erase( get<2>(nodeListOfSols[k]).begin() + removePosition);
                get<0>(nodeListOfSols[k]).second = sequenceTime(get<1>(nodeListOfSols[k]), get<2>(nodeListOfSols[k]), 
                                                                mJobs, mSetupTimes);

                // out << "Solution: ";
                // out << "[ ";
                // for(int l = 0; l < get<1>(nodeListOfSols[k]).size(); l++){
                //     out << get<1>(nodeListOfSols[k])[l] << " ";
                // }
                // out << "]" << endl;
                // out << "[ ";
                // for(int l = 0; l < get<2>(nodeListOfSols[k]).size(); l++){
                //     out << get<2>(nodeListOfSols[k])[l] << " ";
                // }
                // out << "]" << endl;
                // out << "Idle: " << get<0>(nodeListOfSols[k]).first << endl;
                // out << "CTime: " << get<0>(nodeListOfSols[k]).second << endl << endl;
            }

            
            sort(nodeListOfSols.begin(), nodeListOfSols.end(), compBSLWT);

            for (int k = 0; k < nodeListOfSols.size(); k++){
                levelListOfSols.push_back(nodeListOfSols[k]);  
                if(k == totalBranchesPerNode-1) break; 
            } 
        }  

        if(levelListOfSols.size() > totalBranchesPerLevel)
            sort(levelListOfSols.begin(), levelListOfSols.end(), compBSLCT);

        // out << "Nós avaliados: " << levelListOfSols.size() << endl;

        // for(int l = 0; l < levelListOfSols.size(); l++){
        //     out << "[ ";
        //     for(int k = 0; k < get<1>(levelListOfSols[l]).size(); k++){
        //         out << get<1>(levelListOfSols[l])[k] << " ";
        //     }
        //     out << "]" << " Idle: " << get<0>(levelListOfSols[l]).first << 
        //     " CTime: " << get<0>(levelListOfSols[l]).second << endl;        
        // }

        // out << endl;

        listOfSols.clear();
        for (int j = 0; j < levelListOfSols.size(); j++){
            listOfSols.push_back(levelListOfSols[j]);
            if(j == totalBranchesPerLevel) break;
        }

        cout << "Best: " << get<0>(listOfSols[0]).second << endl << endl;
        // cout << "Costs: [ ";
        // for (int j = 0; j < listOfSols.size(); j++){
        //     cout << get<0>(listOfSols[j]) << " ";
        //     if(j == 10) break;
        // }
        // cout << " ]" << endl << endl;
    }

    bestSol.first = get<1>(listOfSols[0]);
    bestSol.second = get<0>(listOfSols[0]).second;
}

void beamSearchDepth(pair< vector<int>, double> &bestSol, pair< vector<int>, double> sol, vector<int> remainigVertices, 
                vector<int> &branchesPerLevel, int totalBranchesPerNode, int totalBranchesPerLevel, int currentLevel){
    // Retornar caso:
    // Limite máximo de branching por nível for atingido ou
    // Não há mais novos vertices a serem adicionados
    // Caso a solução seja melhor que a melhor solução já encontrada, substitua

    // printSolution(remainigVertices, mJobs, mSetupTimes);
    if(remainigVertices.size() == 0){
        if(sol.second < bestSol.second) bestSol = sol;
        return;
    } else if(branchesPerLevel[currentLevel] >= totalBranchesPerLevel){
        return;
    }
    
    vector< pair<vector<int>, double> >  allPossibleInsertions;
    vector<int> newRemainingVertices;
    pair<vector<int>, double> possibleInsertion;
    int inseritonPosition = sol.first.size() - 1;
    int remainingVerticesQtd = remainigVertices.size();
    int removePosition;
    branchesPerLevel[currentLevel]++;
    currentLevel++;

    // Testar todas as possibilidades possíveis para a próxima inserção
    for (int i = 0; i < remainingVerticesQtd; i++){
        allPossibleInsertions.push_back(possibleInsertion);
        allPossibleInsertions[i].first = sol.first;
        allPossibleInsertions[i].first.push_back(remainigVertices[i]);
        // allPossibleInsertions[i].second = sequenceTime(allPossibleInsertions[i].first, mJobs, mSetupTimes).first;
    }
    //cout << "0" << endl;
    // Escolher as N melhores possibilidades de inserção
    sort(allPossibleInsertions.begin(), allPossibleInsertions.end(), compBS);

    //cout << "1" << endl;
    // Realize N branchings, respeitando o espaço por nível
    if(remainigVertices.size() < totalBranchesPerNode) {
        for (size_t i = 0; i < remainigVertices.size(); i++){
            newRemainingVertices = remainigVertices;
            //cout << "1.1" << endl;
            //printSolution(newRemainingVertices, mJobs, mSetupTimes);
            //printSolution(allPossibleInsertions[i].first, mJobs, mSetupTimes);
            removePosition = binarySearch(newRemainingVertices, 0, newRemainingVertices.size()-1,
                                        allPossibleInsertions[i].first[ allPossibleInsertions[i].first.size()-1 ]);
            //cout << "1.2" << endl;
            newRemainingVertices.erase(newRemainingVertices.begin() + removePosition);
            //printSolution(newRemainingVertices, mJobs, mSetupTimes); 
            //cout << "1.3" << endl;
            beamSearchDepth(bestSol, allPossibleInsertions[i], newRemainingVertices, branchesPerLevel, totalBranchesPerNode, 
                    totalBranchesPerLevel, currentLevel); 
            
            //printSolution(newRemainingVertices, mJobs, mSetupTimes); 

            //cout << "1.4" << endl;
        }
    } else {
        for (size_t i = 0; i < totalBranchesPerNode; i++){
            newRemainingVertices = remainigVertices;
            //cout << "1.1" << endl;
            //printSolution(newRemainingVertices, mJobs, mSetupTimes);
            //printSolution(allPossibleInsertions[i].first, mJobs, mSetupTimes);
            removePosition = binarySearch(newRemainingVertices, 0, newRemainingVertices.size()-1,
                                        allPossibleInsertions[i].first[ allPossibleInsertions[i].first.size()-1 ]);
            //cout << "1.2" << endl;
            newRemainingVertices.erase(newRemainingVertices.begin() + removePosition);
            //printSolution(newRemainingVertices, mJobs, mSetupTimes); 
            //cout << "1.3" << endl;
            beamSearchDepth(bestSol, allPossibleInsertions[i], newRemainingVertices, branchesPerLevel, totalBranchesPerNode, 
                    totalBranchesPerLevel, currentLevel); 
            
            //printSolution(newRemainingVertices, mJobs, mSetupTimes); 

            //cout << "1.4" << endl;
        }
    }
    
    //cout << "2" << endl;
    
};

int binarySearch(vector<int> &vec, int posInit, int posEnd, int compare){
    //cout << "Binary Search: " << endl;
    //cout << posInit + posEnd << endl;
    int midlePos = posInit + posEnd;
    midlePos = (midlePos == 0) ? 0 : midlePos / 2;
    //cout << "1.1.1" << endl;
    //cout << vec[midlePos] << endl;
    //cout << endl;
    if(vec[midlePos] == compare){
        return midlePos;
    } else if(posInit == posEnd || posInit > posEnd){
        return -1;
    } else if(compare > vec[midlePos]){
        if(vec[posEnd] == compare) return posEnd;
        binarySearch(vec, midlePos+1, posEnd-1, compare);
    } else {
        if(vec[posInit] == compare) return posInit;
        binarySearch(vec, posInit+1, midlePos+1, compare);
    }
}

bool compBS(pair< vector<int>, double > a, pair< vector<int>, double > b) {
    return a.second < b.second;
}

bool compBSLWT(tuple <pair<double, double>, vector<int>, vector<int> > a, tuple <pair<double, double>, vector<int>, vector<int> > b) {
    return get<0>(a).first*get<0>(a).second < get<0>(b).first*get<0>(b).second;
    //return get<0>(a).first < get<0>(b).first;
}

bool compBSLCT(tuple <pair<double, double>, vector<int>, vector<int> > a, tuple <pair<double, double>, vector<int>, vector<int> > b) {
    return get<0>(a).second < get<0>(b).second;
}

//End

void printData(int n, double ** mJobs, double ** mSetupTimes){
    cout <<"Jobs QTD: "<< n << "\n";

    for (int i = 1; i <= n; i++){
        for(int j = 0; j < 4; j++){
            cout << mJobs[i][j] << " ";
        }
    }
    
   for (int i = 0; i <= n; i++){
        for(int j = 1; j <= n; j++){
            cout << mSetupTimes[i][j] << " ";
        }
    }
 } 

vector<int> construction(int n, double ** mJobs, double ** mSetupTimes, double alfa, double &cost){
    vector<int> s;
    vector<int> listaCandidatos;
    double custoAtual = 0;
    //double exetime;

    for(int i = 1; i <= n; i++){
        listaCandidatos.push_back(i); // insere todos os nós na lista de candidatos
    }


    for(int j = 0; j < 1; j++){ // tamanho subsequencia inicial
        int k = rand() % listaCandidatos.size();
        //s.insert(s.begin() + 1, listaCandidatos[k]); // adiciona trabalho aleatoria a solução
        s.push_back(listaCandidatos[k]);
        listaCandidatos.erase(listaCandidatos.begin() + k); // apaga da lista de candidatos oq ja foi pra solução
        if(j == 0){
            custoAtual += mJobs[s[0]][1] + mSetupTimes[s[0]][s[j]] + mJobs[s[0]][2]; 
        }
        else{
            if(custoAtual >= mJobs[s[j]][1]){
                custoAtual += mSetupTimes[s[j-1]][s[j]] + mJobs[s[j]][2];
            }
            else{
                custoAtual = mSetupTimes[s[j-1]][s[j]] + mJobs[s[j]][2] + mJobs[s[j]][1];
            }
        }
    }

    vector<CustoIn> custoInsertion = calculaCusto(listaCandidatos, s, custoAtual); // calcula os custo de inserção dos candidatos
    std::sort(custoInsertion.begin(), custoInsertion.end(), comp); // ordena de forma crescente de acordo com os custos
    //cout << "print1\n";
    int sel;
    //cout << "tam" << custoInsertion.size() << "\n";
    while(!listaCandidatos.empty()){
        if(alfa == 0){
            sel = 0;
        }
        else{
            sel = rand() % ((int)std::floor(alfa * (custoInsertion.size() - 1)) + 1); // escolhe um nó dentro de uma faixa definida por alfa
        }
        //cout << "print2\n";
        //s.insert(s.begin() + custoInsertion[sel].arestaOut, custoInsertion[sel].noIn); // insere o nó na solucao
        s.push_back(custoInsertion[sel].noIn);
        //cout << custoInsertion[sel].noIn << endl;
        custoAtual += max(custoInsertion[sel].custo,0.0) + mSetupTimes[s[s.size()-2]][s[s.size()-1]] + mJobs[s[s.size()-1]][2]; 
        //cout << "print3\n";
        for(int i = 0; i < listaCandidatos.size(); i++){
            if(listaCandidatos[i]==custoInsertion[sel].noIn)
                listaCandidatos.erase(listaCandidatos.begin() + i); // exclui o nó da lista de candidatos
            }
        
        custoInsertion.erase(custoInsertion.begin(), custoInsertion.end()); // exclui o nó da lista de custos
        custoInsertion = calculaCusto(listaCandidatos, s, custoAtual); // calcula os novos custos de inserção
        std::sort(custoInsertion.begin(), custoInsertion.end(), comp); // ordena os custos
        
    }
    
    setSequencesMatrix(sequencesMatrix,s,n,mJobs,mSetupTimes);
    //compTime = sequenceTime(s, mJobs, mSetupTimes);
    infoSequence initialSolution;
    initialSolution = concatenateSequencev2(mSetupTimes, mJobs, dJob, sequencesMatrix[0][n-1]);
    //idleTime = initialSolution.waitingTime;
    cost = initialSolution.initialTime + initialSolution.duration;
    //printSolution(s, mJobs, mSetupTimes);

    return s;
}

vector<CustoIn> calculaCusto(vector<int> listaCandidatos, vector<int> &s, double custoAtual){
  vector<CustoIn> custoInsertion (listaCandidatos.size());
    //for(int i = 0, l = 0; i < s.size()-1; i++){
        int l = 0;
        for(auto k : listaCandidatos){
            custoInsertion[l].custo = mJobs[k][1] - custoAtual;
            custoInsertion[l].noIn = k; // nó inserido
            //custoInsertion[l].arestaOut = i; // posicao de inserção;
            l++;
        }

  return custoInsertion;

}

double idleTime(vector<int> s, double cTime, int node, double ** mJobs, double **mSetupTimes){
    double totalWT;
    if(s.size() > 0) {
        cTime = mSetupTimes[0][s[0]] + mJobs[s[0]][2] + mJobs[s[0]][1];
        //Calcula o makespan parcial
        for(int i = 0, j = 1; j < s.size(); i++ ,j++){
            // out << i << " " << j << endl;
            if(cTime >= mJobs[s[j]][1]){
                cTime += mSetupTimes[s[i]][s[j]] + mJobs[s[j]][2];
            } else {
                cTime += mSetupTimes[s[i]][s[j]] + mJobs[s[j]][2] + (mJobs[s[j]][1] - cTime);
            }
        }
    } else {
        cTime = 0;
    }

    double aux = mJobs[node][1] - cTime;
    int lastNode = (s.size() > 0) ? s[s.size()-1] : 0;
    double setupTime = mSetupTimes[lastNode][node];
    
    //out << "RD: " << mJobs[node][1] << " - CT: " << cTime << " - ST: " << setupTime << " - LastNode: " << lastNode << endl;
    totalWT = (aux > 0) ? aux + setupTime : setupTime;

    return totalWT;
}

double sequenceTime(vector<int> s, vector<int> remainingVertices, double ** mJobs, double **mSetupTimes){
    double cTime = mSetupTimes[0][s[0]] + mJobs[s[0]][2] + mJobs[s[0]][1];
    double minimunReleaseDate = numeric_limits<double>::max();
    double remainingProcessingTime = 0;
    double remainingLBSetupTime = 0;
    double minimumSetupTime = numeric_limits<double>::max();
    vector<int> unionVertices = remainingVertices;
    unionVertices.push_back(s[s.size()-1]);

    //Calcula o makespan parcial
    for(int i = 0, j = 1; j < s.size(); i++ ,j++){
        // out << i << " " << j << endl;
        if(cTime >= mJobs[s[j]][1]){
            cTime += mSetupTimes[s[i]][s[j]] + mJobs[s[j]][2];
        } else {
            cTime += mSetupTimes[s[i]][s[j]] + mJobs[s[j]][2] + (mJobs[s[j]][1] - cTime);
        }
    }

    //cout << "Makespan: " << cTime << endl;

    //Calcula menor release date
    for(int i = 0; i < remainingVertices.size(); i++){
        if(mJobs[remainingVertices[i]][1] < minimunReleaseDate) minimunReleaseDate = mJobs[remainingVertices[i]][1];
    }
    if( minimunReleaseDate == numeric_limits<double>::max()) minimunReleaseDate = 0;
    //cout << "Minimun RD: " << minimunReleaseDate << endl;

    //Calcula o processing time restante
    for(int i = 0; i < remainingVertices.size(); i++){
        remainingProcessingTime += mJobs[remainingVertices[i]][2];
    }

    //cout << "Remaining PT: " << remainingProcessingTime << endl;

    //Calculando LB setup times
    for(int i = 0; i < remainingVertices.size(); i++){
        for(int j = 0; j < unionVertices.size(); j++){
            if(unionVertices[j] == remainingVertices[i]) continue;
            if(mSetupTimes[unionVertices[j]][remainingVertices[i]] < minimumSetupTime)
                minimumSetupTime = mSetupTimes[unionVertices[j]][remainingVertices[i]];
        }
        remainingLBSetupTime += minimumSetupTime;
        minimumSetupTime = numeric_limits<double>::max();
    }

    //cout << "LB ST: " << remainingLBSetupTime << endl;

    cTime = (cTime > minimunReleaseDate) ? cTime : minimunReleaseDate;
    cTime += remainingProcessingTime;
    cTime += remainingLBSetupTime;

    //cout << "LB Makespan: " << cTime << endl;

    return cTime;
} 

bool comp(const CustoIn& a, const CustoIn& b) // comparação dos custos utilizada para ordenar os objetos
{
    return a.custo < b.custo;
}

double compCostSwap2(vector<int> &s, int i, int j, infoSequence **sequencesMatrix){
 double cost;
    infoSequence swaproute;
    infoSequence dJob;

    dJob.initialTime = 0;
    dJob.duration = 0;
    dJob.firstJob = 0;
    dJob.lastJob = 0;
  
    if(i == 0){
        if(j == 1){
            //dJob.initialTime = 0;
            infoSequence dummyJob = concatenateSequencev2(mSetupTimes, mJobs, dJob, sequencesMatrix[j][i]);
            //swaproute = concatenateSequence(mSetupTimes, mJobs, sequencesMatrix[j][i], sequencesMatrix[j+1][n-1]);
            swaproute = concatenateSequencev2(mSetupTimes, mJobs, dummyJob, sequencesMatrix[j+1][s.size()-1]);
            cost = swaproute.initialTime + swaproute.duration;
        }
        else{
            //dJob.initialTime = 0;
            infoSequence dummyJob = concatenateSequencev2(mSetupTimes, mJobs, dJob, sequencesMatrix[j][j]);
            swaproute = concatenateSequencev2(mSetupTimes, mJobs, dummyJob, sequencesMatrix[1][j-1]);
            //swaproute = concatenateSequencev2(mSetupTimes, mJobs, msequencesMatrix[j][j], sequencesMatrix[i+1][j-1]);
            swaproute = concatenateSequencev2(mSetupTimes, mJobs, swaproute, sequencesMatrix[i][i]);
            swaproute = concatenateSequencev2(mSetupTimes, mJobs, swaproute, sequencesMatrix[j+1][n-1]);
            cost = swaproute.initialTime + swaproute.duration;
        }
    }
    else{
        if(i == n-2){
            //dJob.initialTime = mJobs[s[0]][1];
            infoSequence dummyJob = concatenateSequencev2(mSetupTimes, mJobs, dJob, sequencesMatrix[0][i-1]);
            swaproute = concatenateSequencev2(mSetupTimes, mJobs, dummyJob, sequencesMatrix[j][i]);
            //swaproute = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[0][i-1], sequencesMatrix[j][i]);
            cost = swaproute.initialTime + swaproute.duration;
        }
        else{
            if(j == i+1){
                //dJob.initialTime = mJobs[s[0]][1];
                infoSequence dummyJob = concatenateSequencev2(mSetupTimes, mJobs, dJob, sequencesMatrix[0][i-1]);
                swaproute = concatenateSequencev2(mSetupTimes, mJobs, dummyJob, sequencesMatrix[j][i]);
                //swaproute = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[0][i-1], sequencesMatrix[j][i]);
                swaproute = concatenateSequencev2(mSetupTimes, mJobs, swaproute, sequencesMatrix[j+1][n-1]);
                cost = swaproute.initialTime + swaproute.duration;
            }
            else{
                //dJob.initialTime = mJobs[s[0]][1];
                infoSequence dummyJob = concatenateSequencev2(mSetupTimes, mJobs, dJob, sequencesMatrix[0][i-1]);
                swaproute = concatenateSequencev2(mSetupTimes, mJobs, dummyJob, sequencesMatrix[j][j]);
                //swaproute = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[0][i-1], sequencesMatrix[j][j]);
                swaproute = concatenateSequencev2(mSetupTimes, mJobs, swaproute, sequencesMatrix[i+1][j-1]);
                swaproute = concatenateSequencev2(mSetupTimes, mJobs, swaproute, sequencesMatrix[i][i]);
                swaproute = concatenateSequencev2(mSetupTimes, mJobs, swaproute, sequencesMatrix[j+1][n-1]);
                cost = swaproute.initialTime + swaproute.duration;
            }    
        }
    }
    return cost;
}

void swap(vector<int> &solution, double &cost, infoSequence **sequencesMatrix){
  //double inicioSwap = cpuTime();
  double delta, delta2;
  double menor = cost;
  int pos_i = -1, pos_j = -1; // guarda as posições para realizar a operação
    for(int i = 0; i < solution.size() - 1; i++){ // exclui da operação a primeira e a ultima posição do vetor
        for(int j = i + 1; j < solution.size(); j++){
           delta = compCostSwap2(solution, i, j, sequencesMatrix);
            if(delta < menor){
                menor = delta;
                pos_i = i;
                pos_j = j;
            }
        }
    }

    if(pos_i >= 0){ // realiza a operação
        melhoras++;
        swap(solution[pos_i], solution[pos_j]);
        cost =  menor;
        //setSequencesMatrix(sequencesMatrix,solution,n,mJobs,mSetupTimes);
        melhorasSwap++;

        swap(sequencesMatrix[pos_i][pos_i], sequencesMatrix[pos_j][pos_j]);

        for(int i = 0; i < pos_i; i++){
            for(int j = pos_i; j < n; j++){
                sequencesMatrix[i][j] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[i][j-1], sequencesMatrix[j][j]);
                sequencesMatrix[j][i] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[j][j], sequencesMatrix[j-1][i]);
            }
        }

        for(int j = pos_i + 1; j < n; j++){
            sequencesMatrix[pos_i][j] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[pos_i][j-1], sequencesMatrix[j][j]);
            sequencesMatrix[j][pos_i] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[j][j], sequencesMatrix[j-1][pos_i]);
        }
        
        for(int i = pos_i + 1; i < pos_j; i++){
            for(int j = pos_j; j < n; j++){
                sequencesMatrix[i][j] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[i][j-1], sequencesMatrix[j][j]);
                sequencesMatrix[j][i] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[j][j], sequencesMatrix[j-1][i]);
            }
        }  

         for(int j = pos_j + 1; j < n; j++){
            sequencesMatrix[pos_j][j] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[pos_j][j-1], sequencesMatrix[j][j]);
            sequencesMatrix[j][pos_j] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[j][j], sequencesMatrix[j-1][pos_j]);
        }
    
    }
  //double fimSwap = cpuTime();
 // swapTime += (fimSwap - inicioSwap);   
}
 
double compCostReinsertionv2(int l, vector<int> &s, int i, int j, infoSequence **sequencesMatrix){
    infoSequence swaproute;
    infoSequence dJob;
    infoSequence dummyJob;

    dJob.initialTime = 0;
    dJob.duration = 0;
    dJob.firstJob = 0;
    dJob.lastJob = 0;

    if(i < j){
        if(i == 0){
            dummyJob = concatenateSequencev2(mSetupTimes, mJobs, dJob, sequencesMatrix[i+l][j+l-1]);
            swaproute = concatenateSequencev2(mSetupTimes, mJobs, dummyJob, sequencesMatrix[i][i+l-1]);
            if(j != n-l)
                swaproute = concatenateSequencev2(mSetupTimes, mJobs, swaproute, sequencesMatrix[j+l][n-1]);
        }
        else{
            dummyJob = concatenateSequencev2(mSetupTimes, mJobs, dJob, sequencesMatrix[0][i-1]);
            swaproute = concatenateSequencev2(mSetupTimes, mJobs, dummyJob, sequencesMatrix[i+l][j+l-1]);
            swaproute = concatenateSequencev2(mSetupTimes, mJobs, swaproute, sequencesMatrix[i][i+l-1]);
            if(j != n-l)
                swaproute = concatenateSequencev2(mSetupTimes, mJobs, swaproute, sequencesMatrix[j+l][n-1]);
        }
    }
    else{
        if(j == 0){
            dummyJob = concatenateSequencev2(mSetupTimes, mJobs, dJob, sequencesMatrix[i][i+l-1]);
            swaproute = concatenateSequencev2(mSetupTimes, mJobs, dummyJob, sequencesMatrix[0][i-1]);
            if(i != n-l)
                swaproute = concatenateSequencev2(mSetupTimes, mJobs, swaproute, sequencesMatrix[i+l][n-1]);
        }
        else{
            dummyJob = concatenateSequencev2(mSetupTimes, mJobs, dJob, sequencesMatrix[0][j-1]);
            swaproute = concatenateSequencev2(mSetupTimes, mJobs, dummyJob, sequencesMatrix[i][i-1+l]);
            swaproute = concatenateSequencev2(mSetupTimes, mJobs, swaproute, sequencesMatrix[j][i-1]);
            if(i != n-l)
                swaproute = concatenateSequencev2(mSetupTimes, mJobs, swaproute, sequencesMatrix[i+l][n-1]);
        }
       
    }
    double cost = swaproute.initialTime + swaproute.duration;

    return cost;
}

void reInsertion(int l, vector<int> &solucao, double &custo){ // reinsere um nó em posição diferente
  //double inicioreinsertion = cpuTime();
    double menor = custo;
    double delta, delta2;
    int pos_i = -1, pos_j = -1;
    for(int i = 0; i < solucao.size() - l; i++){ // varre a solução exceto o 0 e o final
        for(int j = 0; j < solucao.size() - l; j++){
            if(i != j){ // exclui a posição inicial do nó
                delta = compCostReinsertionv2(l, solucao, i, j, sequencesMatrix);
                if(delta < menor){
                    menor = delta;
                    pos_i = i;
                    pos_j = j;
                }
            }
        }
    }
    if(pos_i != -1){
        //cout << "i: " << pos_i << " j: " << pos_j << endl;     
        vector<int> subsequence(solucao.begin() + pos_i, solucao.begin() + pos_i + l);
        solucao.erase(solucao.begin() + pos_i, solucao.begin() + pos_i + l);
        solucao.insert(solucao.begin() + pos_j, subsequence.begin(), subsequence.end());
        custo = menor;
        //setSequencesMatrix(sequencesMatrix,solucao,n,mJobs,mSetupTimes);
        melhoras++;
        melhorasReinsert[l-1]++;



////// SEQUENCE ATT         
         if(pos_i > pos_j){
            swap(pos_i, pos_j);
        }

        for(int i = pos_i; i < pos_j + l; i++){
            sequencesMatrix[i][i].firstJob = solucao[i];
            sequencesMatrix[i][i].lastJob = solucao[i];
            sequencesMatrix[i][i].duration = mJobs[solucao[i]][2];
            sequencesMatrix[i][i].initialTime = mJobs[solucao[i]][1];
        }
       
        for(int i = 0; i < pos_i; i++){
            for(int j = pos_i; j < n; j++){
                sequencesMatrix[i][j] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[i][j-1], sequencesMatrix[j][j]);
                sequencesMatrix[j][i] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[j][j], sequencesMatrix[j-1][i]);
            }
        }
        for(int i = pos_i; i < pos_j; i++){
            for(int j = i + 1; j < n; j++){
                sequencesMatrix[i][j] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[i][j-1], sequencesMatrix[j][j]);
                sequencesMatrix[j][i] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[j][j], sequencesMatrix[j-1][i]);
            }
        }
        for(int i = pos_j; i < pos_j+l ; i++){
            for(int j = i + 1; j < n; j++){
                sequencesMatrix[i][j] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[i][j-1], sequencesMatrix[j][j]);
                sequencesMatrix[j][i] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[j][j], sequencesMatrix[j-1][i]);
            }
        }
    }
  //double fimReinsertion =  cpuTime();
  //tempo_reinsertion += (fimReinsertion - inicioreinsertion);
}

vector<int> perturb(vector<int> solucao){
  vector<int> s = solucao;

  int tam_max = ceil(((double)n)/10);
  //int tam_max = 15; // tamanho maximo das subsequencias
  int inicio1, fim1, inicio2, fim2;

  inicio1 = (rand() % ((n - (2 * tam_max)) - 1 + 1)); // posicao minima = 0, posicao maxima = final - 2*tmax
  fim1 = (rand() % ((inicio1 + (tam_max - 1)) - (inicio1 + 1) + 1)) + (inicio1 + 1); // minima = inicio+1, maxima = inicio1 + tmax - 1;
  inicio2 = (rand() % ((n - tam_max) - (fim1 + 1) + 1) + (fim1 + 1)); // minima = fim1 + 1, maxima = final - tmax;
  fim2 = (rand() % ((inicio2 + (tam_max - 1)) - (inicio2 + 1) + 1) + (inicio2 + 1)); // minima = inicio2 + 1, maxima = inicio2 + tmax - 1;


  int d1 = fim1 - inicio1; // tamanho da primeira subsequencia, usado para corrigir erros na inserção
  int d2 = fim2 - inicio2; // tamanho da segunda subsequencia, usado pra corrigir erros na inserção


  //cout << "inicio1: " << inicio1 <<  " fim2: " << fim2 << endl;

  s.erase(s.begin() + inicio2, s.begin() + fim2 + 1); // apaga primeira subsequencia
  s.erase(s.begin() + inicio1, s.begin() + fim1 + 1); // apaga segunda subsequencia
  s.insert(s.begin() + inicio1, &solucao[inicio2], &solucao[fim2] + 1);  // inclui a segunda subsequencia na posicao da primeira
  s.insert(s.begin() + inicio2 + (-1*(d1 - d2)) , &solucao[inicio1], &solucao[fim1] + 1); // inclui a segunda subsequencia na posicao da segunda

  //compCompletionTime(s,compTimes);
  setSequencesMatrix(sequencesMatrix,s,n,mJobs,mSetupTimes);

  return s;
}

void rvnd(vector<int> &solucao, double &custo){
  vector<int> s = solucao;
  //vector<int> nLista = {0,1,2,3,4,5,6,7,8,9,10,11,12,13}; // lista de estruturas in02_360.dat,16876,16902.2,67.8574,72.3509
  //vector<int> nLista = {0,1,2,3}; //instances/in02_360.dat,17081,17101.1,28.9016,32.833
  vector<int> nLista = listaSub;
  double custoMod =  custo;
  int sel, pos;

  while(!nLista.empty()){ // roda enquanto existirem estruturas de vizinhança na lista

    int k = rand() % nLista.size();

    switch(nLista[k]){
        case 0:
            swap(solucao, custoMod, sequencesMatrix);
            break;

        case 1:
            reInsertion(1,solucao, custoMod);
            break;

        case 2:
            reInsertion(2,solucao, custoMod);      
            break;

        case 3:
            reInsertion(3,solucao, custoMod);
            break;
    
        case 4:
            reInsertion(4,solucao, custoMod);   
            break;

        case 5:
            reInsertion(5,solucao, custoMod);        
            break;

        case 6:
            reInsertion(6,solucao, custoMod);       
            break;

        case 7:
            reInsertion(7,solucao, custoMod);  
            break;
    
        case 8:
            reInsertion(8,solucao, custoMod);  
            break;

        case 9:
            reInsertion(9,solucao, custoMod);  
            break;

        case 10:
            reInsertion(10,solucao, custoMod);  
            break;
        
        case 11:
            reInsertion(11,solucao, custoMod);
            break;
        
        case 12:
            reInsertion(12,solucao, custoMod);
            break;
        
        case 13:
            reInsertion(13,solucao, custoMod);
            break;

      

    }

    //custoMod = custoTotal(solucao); // calcula o custo do Movimento

    if(custo > custoMod){ // movimento melhorou o custo
      custo = custoMod;
      s = solucao;
      //nLista = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
      //nLista = {0,1,2,3};
      nLista = listaSub;
    }
    else { // nao melhorou, exclui o movimento da lista
      solucao = s;
      nLista.erase(nLista.begin() + k);
      custoMod = custo;
    }
  }
}

double cpuTime() {
	static struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	return ((double)usage.ru_utime.tv_sec)+(((double)usage.ru_utime.tv_usec)/((double)1000000));
}

void printSolution(vector<int> solucao, double **mJobs, double ** mSetupTimes){
    cout << "Solution: [ ";
    for(int i = 0; i < solucao.size(); i++){
        cout << solucao[i] << " ";
    }
    cout << " ]" << endl;
    /*for(int i = 0; i < n; i++){
        cout << mJobs[solucao[i]][1] << " ";
    }*/

    double cTime = mSetupTimes[0][solucao[0]] + mJobs[solucao[0]][2] + mJobs[solucao[0]][1];
    double totalWT = mJobs[solucao[0]][1];
    /*for(int i = 0; i < s.size(); i++){
        cout << s[i] << " ";
    }*/
    //cout << totalWT;
    for(int i = 0, j = 1; j < solucao.size(); i++ ,j++){
        //cout<< i << " " << j << endl;
        if(cTime >= mJobs[solucao[j]][1]){
            cTime += mSetupTimes[solucao[i]][solucao[j]] + mJobs[solucao[j]][2];
            //cout << " " << 0;
        }
        else{
            totalWT = mJobs[solucao[j]][1] - cTime;
            //cout << " " << totalWT;
            cTime += mSetupTimes[solucao[i]][solucao[j]] + mJobs[solucao[j]][2] + (mJobs[solucao[j]][1] - cTime);
        }
    }

    cout << "TotalWT: " <<  totalWT << endl;
    cout << "CTime: " <<  cTime << endl;
}