/*
 * (mu+lambda)-ES-parallel.h
 *
 *  Created on: 28 lug 2017
 *      Author: RedShy
 */

#ifndef SRC__MU_LAMBDA__ES_THREADS_H_
#define SRC__MU_LAMBDA__ES_THREADS_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <queue>
#include <thread>
#include "ES_MatchingSchema.h"
#include "EditDistance.h"
#include "MatchingSchema.h"

#define CLOCKS_PER_MS (CLOCKS_PER_SEC / 1000)

void evolutionStrategy_t(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, const size_t& p2, matching_schema<bool>& m, edit_distance& e,
		unsigned max_children, ES_MatchingSchema* const parents, const unsigned * const blocksig1, const unsigned * const blocksig2,
		const unsigned worstParentCostValue, const unsigned mu, const unsigned lambda,unsigned offSetInPool, const unsigned threadNumber, const unsigned numberOfThreads);

void evolutionStrategy_t_mu(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, const size_t& p2, matching_schema<bool>& m, edit_distance& e,
		unsigned parentsToComputePerThread, ES_MatchingSchema* const parents, const unsigned * const blocksig1, const unsigned * const blocksig2,
		const unsigned mu, unsigned offSetInPool, ES_MatchingSchema startingMS, const unsigned threadNumber, const unsigned numberOfThreads);

int evolutionStrategy_p(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, const size_t& p2, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations, const unsigned mu,
		const unsigned lambda, unsigned numberOfThreads)
{

	clock_t start = clock();
	long double msElapsed = 0;



	//Initialize stuff for the mutator swap2-E
	const unsigned * const blocksig1 = initializeBlocksSwap2E(sig1, p1);
	const unsigned * const blocksig2 = initializeBlocksSwap2E(sig2, p2);

	ES_MatchingSchema startingMS(sig1, sig2);

	ES_MatchingSchema best;
	best.costValue = std::numeric_limits<unsigned int>::max();

	//max hardware cores
	if(numberOfThreads==0)
	{
		numberOfThreads=std::thread::hardware_concurrency();
	}

	//array of handlers of threads for joining them
	std::thread** threads=new std::thread*[numberOfThreads];

	//Generate mu random individuals
	ES_MatchingSchema* const parents = new ES_MatchingSchema[mu + lambda];

	unsigned offSetInPoolThread=0;
	const unsigned parentsToComputePerThread=mu/numberOfThreads;
	for(unsigned thread=0; thread<numberOfThreads; ++thread)
	{
		std::thread* t=new std::thread(evolutionStrategy_t_mu, s1, s2, s1l, s2l, sig1,
				sig2, sig1l, sig2l, p1, p2, std::ref(m), std::ref(e), parentsToComputePerThread, parents, blocksig1, blocksig2,
				 mu, offSetInPoolThread, startingMS, thread, numberOfThreads);

		threads[thread]=t;

		offSetInPoolThread+=parentsToComputePerThread;
	}

	for(unsigned thread=0; thread<numberOfThreads; ++thread)
	{
		threads[thread]->join();
		delete threads[thread];
	}

	const unsigned last = mu - 1;

	//Select the worst parent in the pool
	unsigned worstParentCostValue=parents[0].costValue;
	for(unsigned i=1; i<mu; ++i)
	{
		if(parents[i].costValue > worstParentCostValue)
		{
			worstParentCostValue=parents[i].costValue;
		}
	}

	//max number of children that a thread can generate.
	//the algorithm don't compute a number of children equal to the remainder of lambda/numberOfThreads
	const unsigned maxChildrenPerThread=lambda/numberOfThreads;

	unsigned generation = 0;
	while (generation <= max_generations)
	{
		//to each thread is given a piece of the genetic pool where store the children
		unsigned offSetInPoolThread=0;
		for(unsigned thread=0; thread<numberOfThreads; ++thread)
		{
			std::thread* t=new std::thread(evolutionStrategy_t, s1, s2, s1l, s2l, sig1,
					sig2, sig1l, sig2l, p1, p2, std::ref(m), std::ref(e), maxChildrenPerThread, parents, blocksig1, blocksig2,
					 worstParentCostValue, mu, lambda,offSetInPoolThread,thread,numberOfThreads);

			threads[thread]=t;

			offSetInPoolThread+=maxChildrenPerThread;
		}

		for(unsigned thread=0; thread<numberOfThreads; ++thread)
		{
			threads[thread]->join();
			delete threads[thread];
		}

		//sorting for selecting the best mu individuals and at the same time get the worst parent
		std::sort(parents, parents+mu+lambda);

		worstParentCostValue = parents[last].costValue;

		//Make a random_shuffle for keeping high entropy
		std::random_shuffle(parents,parents+mu);

		generation++;
	}


	std::sort(parents, parents+mu+lambda);
	best=parents[0];
	std::cout<<best.costValue;

	delete[] blocksig1;
	delete[] blocksig2;

	delete[] parents;
	delete[] threads;

	return best.costValue;
}


void evolutionStrategy_t(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, const size_t& p2, matching_schema<bool>& m, edit_distance& e,
		unsigned max_children, ES_MatchingSchema* const parents, const unsigned * const blocksig1, const unsigned * const blocksig2,
		const unsigned worstParentCostValue, const unsigned mu, const unsigned lambda,unsigned offSetInPool, const unsigned threadNumber, const unsigned numberOfThreads)
{
	const unsigned remainderIndividuals = lambda % numberOfThreads;
	const unsigned reversedPosition=numberOfThreads - threadNumber - 1;
	if(remainderIndividuals != 0)
	{
		//we are distribuiting the remainder of individuals to the threads, so we need to calculate our new offset in the pool and number of parents
		if(reversedPosition<remainderIndividuals)
		{
			max_children++;
			offSetInPool+=remainderIndividuals-1-reversedPosition;
		}
	}

	unsigned addedChildren=0;
	for (unsigned i = 0; i < max_children; i++)
	{
		//Choose random parent
		const unsigned p = rand() % mu;

		//Produce child, in the case parents=1 (like this) just clone
		ES_MatchingSchema child = parents[p];

		//mutate child
		child.swap2_enhanced(blocksig1, blocksig2);

		const int newDistance =
				e.edit_distance_matching_schema_enhanced_with_diagonal(s1,
						s2, s1l, s2l, child.sigma1, child.sigma2, sig1l,
						sig2l, m, worstParentCostValue);

		if (newDistance != -1)
		{
			//The child is better than the worst parent,
			child.costValue = newDistance;

			//so he is added to the pool
			parents[mu+offSetInPool+addedChildren] = child;
			addedChildren++;
		}
	}
}

void evolutionStrategy_t_mu(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, const size_t& p2, matching_schema<bool>& m, edit_distance& e,
		unsigned parentsToComputePerThread, ES_MatchingSchema* const parents, const unsigned * const blocksig1, const unsigned * const blocksig2,
		const unsigned mu, unsigned offSetInPool, ES_MatchingSchema startingMS, const unsigned threadNumber, const unsigned numberOfThreads)
{
	const unsigned remainderIndividuals = mu % numberOfThreads;
	const unsigned reversedPosition=numberOfThreads - threadNumber - 1;

	if(remainderIndividuals != 0)
	{
		//we are distribuiting the remainder of individuals to the threads, so we need to calculate our new offset in the pool and number of parents

		if(reversedPosition<remainderIndividuals)
		{
			parentsToComputePerThread++;
			offSetInPool+=remainderIndividuals-1-reversedPosition;
		}
	}


	for (unsigned i = 0; i < parentsToComputePerThread; ++i)
	{
		startingMS.shuffle();

		startingMS.costValue = e.edit_distance_matching_schema_enhanced(s1, s2,
				s1l, s2l, startingMS.sigma1, startingMS.sigma2, sig1l, sig2l,
				m);

		parents[offSetInPool + i] = startingMS;
	}
}


#endif /* SRC__MU_LAMBDA__ES_THREADS_H_ */
