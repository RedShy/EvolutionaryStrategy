/*
 * (mu+lambda)-ES-shuffle.h
 *
 *  Created on: 05 giu 2017
 *      Author: HantolR
 */

#ifndef SRC__MU_LAMBDA__ES_SHUFFLE_H_
#define SRC__MU_LAMBDA__ES_SHUFFLE_H_

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

//std::vector<ES_MatchingSchema> selectBestIndividuals(const unsigned mu,
//		std::vector<ES_MatchingSchema>& individuals, const bool plusSelection);

int evolutionStrategy_shuffle(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations, const unsigned mu,
		const unsigned lambda)
{
	clock_t start = clock();
	long double msElapsed = 0;

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
			child.swap2();

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

//				//substitute the worst parent in the heap
//				std::pop_heap(parents, parents + mu);
//				parents[last] = child;
//				std::push_heap(parents, parents + mu);

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

//		std::cout<<"PARENTS:\n";
//		for(int i=0; i<mu; i++)
//		{
//			 std::cout<<parents[i].costValue<<" ";
//		}
//		std::cout<<"\n";
//		parents = selectBestIndividuals(mu, parents, plusSelection);

		generation++;
	}

	clock_t timeElapsed = clock() - start;
	msElapsed = timeElapsed / CLOCKS_PER_MS;
	std::cout << msElapsed << " " << best.costValue << "\n";

//	m.print_matching_schema(best.sigma1,best.sigma2);

//	std::cout<<best.costValue;
	return best.costValue;
}

//std::vector<ES_MatchingSchema> selectBestIndividuals(const unsigned mu,
//		std::vector<ES_MatchingSchema>& individuals, const bool plusSelection)
//{
//	//TODO using a heap
//	std::vector<ES_MatchingSchema> bestIndividuals;
//
//	make_heap(individuals.begin(), individuals.end());
//
//	for (unsigned i = 0; i < mu; i++)
//	{
//		pop_heap(individuals.begin(), individuals.end());
//		bestIndividuals.push_back(individuals.back());
//		individuals.pop_back();
//	}
//
//	return bestIndividuals;
//}

#endif /* SRC__MU_LAMBDA__ES_SHUFFLE_H_ */
