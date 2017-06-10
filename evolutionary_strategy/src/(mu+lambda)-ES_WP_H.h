/*
 * (mu+lambda)-ES_WP_H.h
 *
 *  Created on: 09 apr 2017
 *      Author: RedShy
 */

#ifndef SRC__MU_LAMBDA__ES_WP_H_H_
#define SRC__MU_LAMBDA__ES_WP_H_H_


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

int evolutionStrategy_WP_H(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& ms_sig1,
		const std::vector<unsigned>& ms_sig2, const size_t& sig1l,
		const size_t& sig2l,

		const size_t& p1, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations, const unsigned mu,
		const unsigned lambda)
{

	clock_t start = clock();
	long double msElapsed = 0;

	unsigned generation = 0;

	ES_MatchingSchema startingMS(ms_sig1, ms_sig2);
//	startingMS.costValue = e.edit_distance_matching_schema_enhanced(s1, s2, s1l,
//			s2l, startingMS.sigma1, startingMS.sigma2, sig1l, sig2l, m);

//Generate mu individuals
	ES_MatchingSchema parents[mu];
//	parents[0] = startingMS;

	ES_MatchingSchema best;

	for (unsigned i = 0; i < mu; ++i)
	{

//		startingMS.mutate();
		startingMS.shuffle();

		//validate matching schema
		if (ES_isValid(startingMS))
		{
			startingMS.costValue = e.edit_distance_matching_schema_enhanced(s1,
					s2, s1l, s2l, startingMS.sigma1, startingMS.sigma2, sig1l,
					sig2l, m);
			parents[i] = startingMS;

			if (parents[i].costValue < best.costValue)
			{
				best = parents[i];

				clock_t timeElapsed = clock() - start;
				msElapsed = timeElapsed / CLOCKS_PER_MS;
				std::cout << msElapsed << " " << best.costValue << "\n";

			}
		}
		else
		{
			//TODO: not valid, maybe mutate until is valid?
			//repeat iteration
			i--;
		}
	}

	const unsigned last = mu - 1;

	std::make_heap(parents, parents + mu);

//	unsigned same = 0;
	while (generation <= max_generations)
	{

		//Generate lambda children. Only mutation, no recombination
		for (unsigned i = 0; i < lambda; i++)
		{
			//Choose random parent
			unsigned p = rand() % mu;

			//Produce child, in the case parents=1 (like this) just clone
			ES_MatchingSchema child = parents[p];

			//mutate child
			child.swap2();

			//validate child
			if (ES_isValid(child))
			{
				//select the worst parent
				unsigned worstParentCostValue = parents[0].costValue;

				int newDistance =
						e.edit_distance_matching_schema_enhanced_with_diagonal(
								s1, s2, s1l, s2l, child.sigma1, child.sigma2,
								sig1l, sig2l, m, worstParentCostValue);

				if (newDistance != -1)
				{
					//The child is better than the worst parent, so he become a new parent
					child.costValue = newDistance;

					std::pop_heap(parents, parents + mu);
					parents[last] = child;
					std::push_heap(parents, parents + mu);

					if (child.costValue < best.costValue)
					{
						best = child;

						clock_t timeElapsed = clock() - start;
						msElapsed = timeElapsed / CLOCKS_PER_MS;
						std::cout << msElapsed << " " << best.costValue << "\n";
					}
				}

//				else child discarded
			}
			else
			{
				//TODO: not valid, maybe mutate until is valid?
				//repeat iteration
				i--;
			}
		}
		generation++;

	}

	clock_t timeElapsed = clock() - start;
	msElapsed = timeElapsed / CLOCKS_PER_MS;
	std::cout << msElapsed << " " << best.costValue << "\n";


	return best.costValue;

}


#endif /* SRC__MU_LAMBDA__ES_WP_H_H_ */
