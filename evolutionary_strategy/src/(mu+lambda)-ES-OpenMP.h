/*
 * (mu+lambda)-ES-OpenMP.h
 *
 *  Created on: 30 lug 2017
 *      Author: RedShy
 */

#ifndef SRC__MU_LAMBDA__ES_OPENMP_H_
#define SRC__MU_LAMBDA__ES_OPENMP_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <queue>
#include <omp.h>
#include "ES_MatchingSchema.h"
#include "EditDistance.h"
#include "MatchingSchema.h"

#define CLOCKS_PER_MS (CLOCKS_PER_SEC / 1000)

void computeChildren(const unsigned lambda, const unsigned mu,
		const unsigned * const blocksig1, const unsigned * const blocksig2,
		const std::vector<unsigned>& s1, const std::vector<unsigned>& s2,
		const size_t& s1l, const size_t& s2l, const size_t& sig1l,
		const size_t& sig2l, const unsigned worstParentCostValue,
		ES_MatchingSchema* const parents, edit_distance& e,
		matching_schema<bool>& m, const unsigned offSetInPoolThread);

int evolutionStrategy_omp(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, const size_t& p2, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations, const unsigned mu,
		const unsigned lambda)
{
//	clock_t start = clock();
//	long double msElapsed = 0;

//	const unsigned numberOfThreads=std::thread::hardware_concurrency;
	const unsigned numberOfThreads=3;

	//Initialize stuff for the mutator swap2-E
	const unsigned * const blocksig1 = initializeBlocksSwap2E(sig1, p1);
	const unsigned * const blocksig2 = initializeBlocksSwap2E(sig2, p2);

	ES_MatchingSchema startingMS(sig1, sig2);

	ES_MatchingSchema best;
	best.costValue = std::numeric_limits<unsigned int>::max();

	//Generate mu random individuals
	ES_MatchingSchema* const parents = new ES_MatchingSchema[mu + lambda];
	for (unsigned i = 0; i < mu; ++i)
	{
		startingMS.shuffle();

		startingMS.costValue = e.edit_distance_matching_schema_enhanced(s1, s2,
				s1l, s2l, startingMS.sigma1, startingMS.sigma2, sig1l, sig2l,
				m);

		parents[i] = startingMS;

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

	unsigned generation = 0;
	while (generation <= max_generations)
	{
		unsigned childrenInPool=0;


		#pragma omp parallel num_threads(numberOfThreads)
		{
			computeChildren(lambda/numberOfThreads, mu, blocksig1,
					blocksig2, s1, s2, s1l, s2l, sig1l, sig2l, worstParentCostValue,
					parents, e, m, (lambda/numberOfThreads)*omp_get_thread_num());

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



	return best.costValue;
}

void computeChildren(const unsigned lambda, const unsigned mu,
		const unsigned * const blocksig1, const unsigned * const blocksig2,
		const std::vector<unsigned>& s1, const std::vector<unsigned>& s2,
		const size_t& s1l, const size_t& s2l, const size_t& sig1l,
		const size_t& sig2l, const unsigned worstParentCostValue,
		ES_MatchingSchema* const parents, edit_distance& e,
		matching_schema<bool>& m, const unsigned offSetInPoolThread)
{
	unsigned childrenInPool=0;
	//Generate lambda children. Only mutation, no recombination
	for (unsigned i = 0; i < lambda; i++)
	{
		//Choose random parent
		const unsigned p = rand() % mu;
		//Produce child, in the case parents=1 (like this) just clone
		ES_MatchingSchema child = parents[p];
		//mutate child
		child.swap2_enhanced(blocksig1, blocksig2);
		const int newDistance =
				e.edit_distance_matching_schema_enhanced_with_diagonal(s1, s2,
						s1l, s2l, child.sigma1, child.sigma2, sig1l, sig2l, m,
						worstParentCostValue);
		if (newDistance != -1)
		{
			//The child is better than the worst parent,
			child.costValue = newDistance;
			//so he is added to the pool

			parents[mu + offSetInPoolThread + childrenInPool] = child;
			childrenInPool++;
		}
	}
}


#endif /* SRC__MU_LAMBDA__ES_OPENMP_H_ */
