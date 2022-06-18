#include <time.h>
#include <math.h>
#include <string>
#include <iostream>
#include <thread>
#include "genetic-alg.h"

GeneticAlgorithm::GaSolution::GaSolution() {
	// empty
}

void GeneticAlgorithm::GaSolution::InitGeneString(unsigned __int16 n, bool randomise) {

	genes = new unsigned __int32[n];

	if (randomise) {

		for (unsigned __int16 i = 0; i < n; i++) {
			genes[i] = GenRandomNumber();
		}
	}
}

void GeneticAlgorithm::Initialise(__int32 inputPopSize, __int32 inputNumberOfGenerations, __int32 numberOfGenes, double(*fitnessFunction)(unsigned __int32* inGenes), __int32 inputGroupSize, __int32 mutationRateIn100000, __int32 crossoverRateIn100000) {
	popSize = inputPopSize;
	generations = inputNumberOfGenerations;
	nGenes = numberOfGenes;
	FitnessEval = fitnessFunction;
	groupSize = inputGroupSize;
	muRate = mutationRateIn100000;
	crRate = crossoverRateIn100000;

	block = false;
	srand((unsigned)time(NULL));
}

void GeneticAlgorithm::CreatePopulation() {
	population = new GaSolution[popSize];

	for (unsigned __int32 i = 0; i < popSize; i++) {
		population[i].InitGeneString(nGenes, true);
	}

	bestSolution.InitGeneString(nGenes, true);
	bestSolution.fitness = FitnessEval(bestSolution.genes);
}

void GeneticAlgorithm::FeThread(unsigned __int32 start, unsigned __int32 end) {

	// race condition possible where UpdateBest is called simultaneously,
	// but the consequence of it occuring is basically nothing and the chance of occurrance is low, so it's not an issue

	for (unsigned __int32 i = start; i < end; i++) {
		population[i].fitness = FitnessEval(population[i].genes);

		if (population[i].fitness > bestSolution.fitness) {
			UpdateBest(population[i]);
		}
	}
}

void GeneticAlgorithm::EvaluateFitnessForPop() {
	if (popSize >= 2000) {
		std::thread t1(FeThread, 0, popSize / 4);
		std::thread t2(FeThread, popSize / 4, popSize / 2);
		std::thread t3(FeThread, popSize / 2, (3 * popSize) / 4);
		std::thread t4(FeThread, (3 * popSize) / 4, popSize);
		t1.join();
		t2.join();
		t3.join();
		t4.join();
	}
	else {
		std::thread t1(FeThread, 0, popSize);
		t1.join();
	}
}

void GeneticAlgorithm::UpdateBest(GeneticAlgorithm::GaSolution input) {
	bestSolution.fitness = input.fitness;
	for (int i = 0; i < nGenes; i++) {
		bestSolution.genes[i] = input.genes[i];
	}
}

void GeneticAlgorithm::GenerateChild(unsigned __int32 posM, unsigned __int32 posD, GeneticAlgorithm::GaSolution newSol) {

	unsigned __int32 tempRand;

	unsigned __int32 crMask = 0;
	unsigned __int32 crIndex = nGenes + 1;

	unsigned __int32 mMask = 0;
	unsigned __int32 mIndex = nGenes + 1;

	if (GenRandomNumber() % 100000 < crRate) {
		//crossover
		tempRand = 1 + (GenRandomNumber() % (32 * nGenes - 1)); // N between 1 and (32*nGenes - 1)

		crIndex = tempRand / 32;
		crMask = pow(2, (tempRand % 32) ) - 1;
	}

	if (GenRandomNumber() % 100000 < muRate) {
		//mutation
		tempRand = GenRandomNumber() % (32 * nGenes); // N between 0 and (32*nGenes - 1)

		mIndex = tempRand / 32;
		mMask = pow(2, (tempRand % 32) );
	}

	for (unsigned __int32 i = 0; i < nGenes; i++) {
		if (i < crIndex)
			newSol.genes[i] = population[posM].genes[i];
		else if (i > crIndex)
			newSol.genes[i] = population[posD].genes[i];
		else
			newSol.genes[i] = (population[posM].genes[i] & crMask) + (population[posD].genes[i] & ~crMask);

		if (i == mIndex)
			newSol.genes[i] = (newSol.genes[i] & ~mMask) + (~(newSol.genes[i]) & mMask);
	}
}

void GeneticAlgorithm::TsThread(GaSolution* nextGen, unsigned __int32 start, unsigned __int32 end) {
	bool unique;
	unsigned __int32 bestContender;
	unsigned __int32* contenders = new unsigned __int32[groupSize];


	for (unsigned __int32 popMember = start; popMember < end; popMember++) {

		// for m
		for (int i = 0; i < groupSize; i++) {
			contenders[i] = 0;
		}

		for (int i = 0; i < groupSize; i++) {
			unique = false;
			while (!unique) {
				contenders[i] = unsigned __int32(GenRandomNumber() % popSize);
				unique = true;
				for (int j = 0; j < groupSize; j++) {
					if (contenders[i] == contenders[j] && i != j) {
						unique = false;
					}
				}
			}
		}
		//initialised contenders, now find the best 1
		bestContender = contenders[0];
		for (int i = 1; i < groupSize; i++) {
			if (population[contenders[i]].fitness > population[bestContender].fitness) {
				bestContender = contenders[i];
			}
		}

		unsigned __int32 mother = bestContender;

		// for d
		for (int i = 0; i < groupSize; i++) {
			contenders[i] = 0;
		}

		for (int i = 0; i < groupSize; i++) {
			unique = false;
			while (!unique) {
				contenders[i] = unsigned __int32(GenRandomNumber() % popSize);
				unique = !(contenders[i] == mother);
				for (int j = 0; j < groupSize; j++) {
					if (contenders[i] == contenders[j] && i != j) {
						unique = false;
					}
				}
			}
		}

		//initialised contenders, now find the best 1
		bestContender = contenders[0];
		for (int i = 1; i < groupSize; i++) {
			if (population[contenders[i]].fitness > population[bestContender].fitness) {
				bestContender = contenders[i];
			}
		}

		unsigned __int32 father = bestContender;

		GenerateChild(mother, father, nextGen[popMember]);
	}
}

void GeneticAlgorithm::TournamentSelection() {
	GaSolution* nextGen = new GaSolution[popSize];
	for (unsigned __int32 i = 0; i < popSize; i++) {
		nextGen[i].InitGeneString(nGenes, false);
	}

	if (popSize >= 2000) {
		std::thread t1(TsThread, nextGen, 0, popSize / 4);
		std::thread t2(TsThread, nextGen, popSize / 4, popSize / 2);
		std::thread t3(TsThread, nextGen, popSize / 2, (3 * popSize) / 4);
		std::thread t4(TsThread, nextGen, (3 * popSize) / 4, popSize);
		t1.join();
		t2.join();
		t3.join();
		t4.join();
	}
	else {
		std::thread t1(TsThread, nextGen, 0, popSize);
		t1.join();
	}



	for (unsigned __int32 i = 0; i < popSize; i++) {
		for (unsigned __int32 j = 0; j < nGenes; j++) {
			population[i].genes[j] = nextGen[i].genes[j];
		}
		delete[] nextGen[i].genes;
	}
	delete[] nextGen;
}

void GeneticAlgorithm::RunGeneticAlgorithm(bool printOutput) {
	if (block) {
		std::cout << "\nGenAlg instance not initialised properly, use the parametised constructor only!\n";
	}
	CreatePopulation();

	std::cout.precision(15);
	for (int i = 0; i < generations; i++) {
		EvaluateFitnessForPop();
		if (printOutput) std::cout << std::scientific << "Generation " << i + 1 << " completed. Best fitness = " << bestSolution.fitness << std::endl;
		TournamentSelection();
	}
}

unsigned __int32 GeneticAlgorithm::GenRandomNumber() {
	// existing randomisation functions weren't being random enough, not even <random>
	// so we'll literally just randomise bits

	unsigned __int32 n = 0;
	for (int j = 0; j < 4; j++) {
		unsigned __int32 val = rand() % 256;
		n += (val << (j * 8));
	}

	return n;
}

void GeneticAlgorithm::ClearObject() {
	for (unsigned __int32 i = 0; i < popSize; i++) {
		delete[] population[i].genes;
	}
	delete[] population;
	block = true;
}


unsigned __int32 GeneticAlgorithm::popSize;
unsigned __int16 GeneticAlgorithm::generations;
unsigned __int16 GeneticAlgorithm::nGenes;
double(*GeneticAlgorithm::FitnessEval)(unsigned __int32* inputGenes);
unsigned __int16 GeneticAlgorithm::groupSize;
unsigned __int32 GeneticAlgorithm::muRate;
unsigned __int32 GeneticAlgorithm::crRate;

bool GeneticAlgorithm::block = true;

GeneticAlgorithm::GaSolution* GeneticAlgorithm::population;
GeneticAlgorithm::GaSolution GeneticAlgorithm::bestSolution;