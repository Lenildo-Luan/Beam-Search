#ifndef infoSeq_H
#define infoSeq_H

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

using namespace std;

struct infoSequence{
    double duration = 0;    //minimum duration
    double initialTime = 0; // earliest time 
    double waitingTime = 0;
    int firstJob = 0;		// first job of the subsequence
	int lastJob = 0;		// last job of the subsequence
};

infoSequence **sequencesMatrix;
infoSequence dJob;

void setSequencesMatrix(infoSequence **sequencesMatrix, int n, double **mJobs, double **setupTimes);
static infoSequence concatenateSequencev2(double **mSetupTimes, double ** mJobs, infoSequence seq1, infoSequence seq2);

inline static infoSequence concatenateSequencev2(double **mSetupTimes, double ** mJobs, infoSequence seq1, infoSequence seq2){

    double waitingTime;
    infoSequence newSequence;


    if(seq1.initialTime + seq1.duration >= seq2.initialTime)
        waitingTime = 0;
    else{
        waitingTime = 0;
        if(seq1.initialTime + seq1.duration < mJobs[seq2.firstJob][1])
            waitingTime += mJobs[seq2.firstJob][1] - (seq1.initialTime + seq1.duration);
        if(seq1.initialTime + seq1.duration + waitingTime + mSetupTimes[seq1.lastJob][seq2.firstJob]  < seq2.initialTime)
            waitingTime += seq2.initialTime - (seq1.initialTime + seq1.duration + mSetupTimes[seq1.lastJob][seq2.firstJob] + waitingTime);
    }
                  
    newSequence.firstJob = seq1.firstJob;
    newSequence.lastJob = seq2.lastJob;
    newSequence.duration = seq1.duration + seq2.duration + mSetupTimes[seq1.lastJob][seq2.firstJob];
    newSequence.initialTime = seq1.initialTime + waitingTime;
    newSequence.waitingTime = waitingTime;    

    return newSequence;
}

void setSequencesMatrix(infoSequence **sequencesMatrix, vector<int> s, int n, double **mJobs, double **mSetupTimes){

    dJob.initialTime = 0;
    dJob.duration = 0;
    dJob.firstJob = 0;
    dJob.lastJob = 0;


    for(int i = 0; i < n; i++){
        sequencesMatrix[i][i].firstJob = s[i];
        sequencesMatrix[i][i].lastJob = s[i];
        sequencesMatrix[i][i].duration = mJobs[s[i]][2];
        sequencesMatrix[i][i].initialTime = mJobs[s[i]][1];
    }

    

    for(int i = 0; i < n; i++){ // original path
        for(int j = i+1; j < n; j++){
            sequencesMatrix[i][j] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[i][j-1], sequencesMatrix[j][j]);
        }
    }

    for(int i = n-1; i >= 1; i--){ // inverted subsequences
        for(int j = i-1; j >= 0; j--){
            sequencesMatrix[i][j] = concatenateSequencev2(mSetupTimes, mJobs, sequencesMatrix[i][j+1], sequencesMatrix[j][j]);
        }

    }
}

#endif