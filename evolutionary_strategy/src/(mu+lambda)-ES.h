/*
 * (mu+lambda)-ES.h
 *
 *  Created on: 05 giu 2017
 *      Author: HantolR
 */

#ifndef SRC__MU_LAMBDA__ES_H_
#define SRC__MU_LAMBDA__ES_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <queue>
#include "ES_MatchingSchema.h"
#include "EditDistance.h"
#include "MatchingSchema.h"

#define CLOCKS_PER_MS (CLOCKS_PER_SEC / 1000)

int evolutionStrategy(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, const size_t& p2, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations, const unsigned mu,
		const unsigned lambda)
{
	clock_t start = clock();
	long double msElapsed = 0;

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

		if (parents[i].costValue < best.costValue)
		{
			best = parents[i];

			clock_t timeElapsed = clock() - start;
			msElapsed = timeElapsed / CLOCKS_PER_MS;
			std::cout << msElapsed << " " << best.costValue << "\n";
		}

	}

	const unsigned last = mu - 1;

	std::sort(parents, parents + mu);
	unsigned worstParentCostValue = parents[last].costValue;
	std::random_shuffle(parents,parents+mu);

	unsigned generation = 0;
	while (generation <= max_generations)
	{
		unsigned childrenInPool = 0;

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
					e.edit_distance_matching_schema_enhanced_with_diagonal(s1,
							s2, s1l, s2l, child.sigma1, child.sigma2, sig1l,
							sig2l, m, worstParentCostValue);

			if (newDistance != -1)
			{
				//The child is better than the worst parent,
				child.costValue = newDistance;

				//so he is added to the pool
				parents[mu + childrenInPool] = child;
				childrenInPool++;

				if (child.costValue < best.costValue)
				{
					best = child;

					clock_t timeElapsed = clock() - start;
					msElapsed = timeElapsed / CLOCKS_PER_MS;
					std::cout << msElapsed << " " << best.costValue << "\n";
				}
			}
		}

		std::sort(parents, parents+mu+childrenInPool);

		worstParentCostValue = parents[last].costValue;

		std::random_shuffle(parents,parents+mu);

		generation++;
	}

	delete[] blocksig1;
	delete[] blocksig2;

	delete[] parents;

	clock_t timeElapsed = clock() - start;
	msElapsed = timeElapsed / CLOCKS_PER_MS;
	std::cout << msElapsed << " " << best.costValue << "\n";

	m.print_matching_schema(best.sigma1,best.sigma2);

//	std::cout<<best.costValue;
//	std::cout<<"CHIAMATE A EDIT DISTANCE= "<<edit_distance::tentativi<<"\n";
	return best.costValue;
}


#endif /* SRC__MU_LAMBDA__ES_H_ */
