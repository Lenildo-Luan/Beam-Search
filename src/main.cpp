#include "readData.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <tuple>

using namespace std;

double ** mJobs; // matriz de adjacencia
double ** mSetupTimes; // matriz reorganizada;
int n; // quantidade total de vertices

int binarySearch(vector<int> &vec, int posInit, int posEnd, int compare);
void beamSearch(pair< vector<int>, double> &bestSol, int totalBranchesPerNode, int totalBranchesPerLevel);
void printSolution(vector<int> solucao, double **mJobs, double ** mSetupTimes);
bool compBSLWT(tuple <pair<double, double>, vector<int>, vector<int> > a, tuple <pair<double, double>, vector<int>, vector<int> > b);
bool compBSLCT(tuple <pair<double, double>, vector<int>, vector<int> > a, tuple <pair<double, double>, vector<int>, vector<int> > b);
double sequenceTime(vector<int> s, vector<int> remainingVertices, double ** mJobs, double **mSetupTimes);
double idleTime(vector<int> s, double cTime, int node, double ** mJobs, double **mSetupTimes);

int main(int argc, char** argv) {
    //Lê os dados da instância
    readData(argc, argv, &n, &mJobs, &mSetupTimes);

    //Declara variáveis
    pair< vector<int>, double> bestSol, sol;
    bestSol.second = numeric_limits<double>::max();
    sol.first = {};
    sol.second = 0;

    //Chama Beam Search e mostra solução
    beamSearch(bestSol, 5, 500);
    printSolution(bestSol.first, mJobs, mSetupTimes);

    return 0;
}

//Algorítimos
void beamSearch(pair< vector<int>, double> &bestSol, int totalBranchesPerNode, int totalBranchesPerLevel){
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
        //Limpa a lista de nós por nível para nova análise
        levelListOfSols.clear();

        for (int j = 0; j < listOfSols.size(); j++){
            //Limpa a lista de nós para uma nova análise
            nodeListOfSols.clear();
            
            //Obtem todas as possíveis expações a partir do nó J 
            for (int k = 0; k < get<2>(listOfSols[j]).size(); k++){
                nodeListOfSols.push_back(sol); 
                nodeListOfSols[k] = listOfSols[j];
                get<0>(nodeListOfSols[k]).first = idleTime(get<1>(nodeListOfSols[k]), get<0>(nodeListOfSols[k]).second, 
                                                           get<2>(listOfSols[j])[k], mJobs, mSetupTimes);
                get<1>(nodeListOfSols[k]).push_back(get<2>(listOfSols[j])[k]);
                removePosition = binarySearch(get<2>(nodeListOfSols[k]), 0, get<2>(nodeListOfSols[k]).size()-1,
                                              get<1>(nodeListOfSols[k])[ get<1>(nodeListOfSols[k]).size()-1 ]);
                get<2>(nodeListOfSols[k]).erase( get<2>(nodeListOfSols[k]).begin() + removePosition);
                get<0>(nodeListOfSols[k]).second = sequenceTime(get<1>(nodeListOfSols[k]), get<2>(nodeListOfSols[k]), 
                                                                mJobs, mSetupTimes);
            }

            //Ordena pelo Waiting Time a lista de nós obtidos
            sort(nodeListOfSols.begin(), nodeListOfSols.end(), compBSLWT);

            //Escolhe os W melhores nós para expadir
            for (int k = 0; k < nodeListOfSols.size(); k++){
                levelListOfSols.push_back(nodeListOfSols[k]);  
                if(k == totalBranchesPerNode-1) break; 
            } 
        }  

        //Caso a quantidade de nós para expansão ultrapasse o limite por nível, ordena pelo LBMakespan
        if(levelListOfSols.size() > totalBranchesPerLevel)
            sort(levelListOfSols.begin(), levelListOfSols.end(), compBSLCT);

        //Limpa a lista de nós para expanção para nova análise  
        listOfSols.clear();

        //Expande os N melhores nós do nível
        for (int j = 0; j < levelListOfSols.size(); j++){
            listOfSols.push_back(levelListOfSols[j]);
            if(j == totalBranchesPerLevel) break;
        }

    }

    //Obtem melhor solução ao final do algorítimo
    bestSol.first = get<1>(listOfSols[0]);
    bestSol.second = get<0>(listOfSols[0]).second;
}

int binarySearch(vector<int> &vec, int posInit, int posEnd, int compare){
    int midlePos = posInit + posEnd;
    midlePos = (midlePos == 0) ? 0 : midlePos / 2;

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

//Funções de ordenação
bool compBSLWT(tuple <pair<double, double>, vector<int>, vector<int> > a, tuple <pair<double, double>, vector<int>, vector<int> > b) {
    //Ordena pelo Waiting Time

    //Esta comparação apresenta melhoras em algumas instâncias
    //return get<0>(a).first*get<0>(a).second < get<0>(b).first*get<0>(b).second; 
    //Comparação do artigo
    return get<0>(a).first < get<0>(b).first;
}

bool compBSLCT(tuple <pair<double, double>, vector<int>, vector<int> > a, tuple <pair<double, double>, vector<int>, vector<int> > b) {
    //Ordena pelo LB Completion Time
    return get<0>(a).second < get<0>(b).second;
}

//Funçoes de cálculos
double idleTime(vector<int> s, double cTime, int node, double ** mJobs, double **mSetupTimes){
    double idleTime;
    int lastNode = (s.size() > 0) ? s[s.size()-1] : 0;
    double setupTime = mSetupTimes[lastNode][node];

    //Calcula o makespan parcial
    if(s.size() > 0) {
        cTime = mSetupTimes[0][s[0]] + mJobs[s[0]][2] + mJobs[s[0]][1];
        for(int i = 0, j = 1; j < s.size(); i++ ,j++){
            if(cTime >= mJobs[s[j]][1]){
                cTime += mSetupTimes[s[i]][s[j]] + mJobs[s[j]][2];
            } else {
                cTime += mSetupTimes[s[i]][s[j]] + mJobs[s[j]][2] + (mJobs[s[j]][1] - cTime);
            }
        }
    } else {
        cTime = 0;
    }

    cTime = mJobs[node][1] - cTime;

    //Calcula idleTime
    idleTime = (cTime > 0) ? cTime + setupTime : setupTime;

    return idleTime;
}

double sequenceTime(vector<int> s, vector<int> remainingVertices, double ** mJobs, double **mSetupTimes){
    double LBMakespan;
    double cTime = mSetupTimes[0][s[0]] + mJobs[s[0]][2] + mJobs[s[0]][1];
    double minimunReleaseDate = numeric_limits<double>::max();
    double remainingProcessingTime = 0;
    double remainingLBSetupTime = 0;
    double minimumSetupTime = numeric_limits<double>::max();
    vector<int> unionVertices = remainingVertices;
    unionVertices.push_back(s[s.size()-1]);

    //Calcula o makespan parcial
    for(int i = 0, j = 1; j < s.size(); i++ ,j++){
        if(cTime >= mJobs[s[j]][1]){
            cTime += mSetupTimes[s[i]][s[j]] + mJobs[s[j]][2];
        } else {
            cTime += mSetupTimes[s[i]][s[j]] + mJobs[s[j]][2] + (mJobs[s[j]][1] - cTime);
        }
    }

    //Calcula menor release date dos nós ainda não inseridos
    for(int i = 0; i < remainingVertices.size(); i++){
        if(mJobs[remainingVertices[i]][1] < minimunReleaseDate) minimunReleaseDate = mJobs[remainingVertices[i]][1];
    }
    if( minimunReleaseDate == numeric_limits<double>::max()) minimunReleaseDate = 0;

    //Calcula o processing time total restante
    for(int i = 0; i < remainingVertices.size(); i++){
        remainingProcessingTime += mJobs[remainingVertices[i]][2];
    }

    //Calcula o lower bound dos setup times
    for(int i = 0; i < remainingVertices.size(); i++){
        for(int j = 0; j < unionVertices.size(); j++){
            if(unionVertices[j] == remainingVertices[i]) continue;
            if(mSetupTimes[unionVertices[j]][remainingVertices[i]] < minimumSetupTime)
                minimumSetupTime = mSetupTimes[unionVertices[j]][remainingVertices[i]];
        }
        remainingLBSetupTime += minimumSetupTime;
        minimumSetupTime = numeric_limits<double>::max();
    }

    //Calcula o lower bound do makespan 
    LBMakespan = (cTime > minimunReleaseDate) ? cTime : minimunReleaseDate;
    LBMakespan += remainingProcessingTime;
    LBMakespan += remainingLBSetupTime;

    return LBMakespan;
} 

//Funções de utilidades
void printSolution(vector<int> solucao, double **mJobs, double ** mSetupTimes){
    cout << "Solution: [ ";
    for(int i = 0; i < solucao.size(); i++){
        cout << solucao[i] << " ";
    }
    cout << " ]" << endl;

    double cTime = mSetupTimes[0][solucao[0]] + mJobs[solucao[0]][2] + mJobs[solucao[0]][1];
    double totalWT = mJobs[solucao[0]][1];

    cout << totalWT;
    for(int i = 0, j = 1; j < solucao.size(); i++ ,j++){
        if(cTime >= mJobs[solucao[j]][1]){
            cTime += mSetupTimes[solucao[i]][solucao[j]] + mJobs[solucao[j]][2];
            cout << " " << 0;
        }
        else{
            totalWT = mJobs[solucao[j]][1] - cTime;
            cout << " " << totalWT;
            cTime += mSetupTimes[solucao[i]][solucao[j]] + mJobs[solucao[j]][2] + (mJobs[solucao[j]][1] - cTime);
        }
    }

    cout << endl << "TotalWT: " <<  totalWT << endl;
    cout << "CTime: " <<  cTime << endl;
}
