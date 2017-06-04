/*
 * (mu+lambda)-ES_WP.h
 *
 *  Created on: 04 apr 2017
 *      Author: RedShy
 */

#ifndef mu_lambda_ES_WP
#define mu_lambda_ES_WP

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

int evolutionStrategy_WP(const std::vector<unsigned>& s1,
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

	//Keep the best only for printing intermediate results, is useless otherwise
	ES_MatchingSchema best;
	best.costValue=std::numeric_limits<unsigned int>::max();

	//Generate mu random individuals
	// TODO: not compatible with g++ version < 5
	//ES_MatchingSchema parents[mu];
	std::vector<ES_MatchingSchema> parents(mu);
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
	//std::make_heap(parents, parents + mu);
	std::make_heap(parents.begin(), parents.end());

	unsigned generation = 0;
	while (generation <= max_generations)
	{
		//Generate lambda children. Only mutation, no recombination
		for (unsigned i = 0; i < lambda; i++)
		{
			//Choose random parent
			unsigned p = rand() % mu;

			//Produce child: just clone the parent
			ES_MatchingSchema child = parents[p];

			//mutate child
			child.swap2();

			//select the worst parent: is the first element of the heap
			unsigned worstParentCostValue = parents[0].costValue;

			int newDistance =
					e.edit_distance_matching_schema_enhanced_with_diagonal(s1,
							s2, s1l, s2l, child.sigma1, child.sigma2, sig1l,
							sig2l, m, worstParentCostValue);

			if (newDistance != -1)
			{
				//The child is better than the worst parent, so he become a new parent
				child.costValue = newDistance;

				//substitute the worst parent in the heap
				//std::pop_heap(parents, parents + mu);
				std::pop_heap(parents.begin(), parents.end());
				parents[last] = child;
				//std::push_heap(parents, parents + mu);
				std::push_heap(parents.begin(), parents.end());

				if (child.costValue < best.costValue)
				{
					best = child;

					clock_t timeElapsed = clock() - start;
					msElapsed = timeElapsed / CLOCKS_PER_MS;
					std::cout << msElapsed << " " << best.costValue << "\n";
				}
			}
		}
		generation++;
	}

	clock_t timeElapsed = clock() - start;
	msElapsed = timeElapsed / CLOCKS_PER_MS;
	std::cout << msElapsed << " " << best.costValue << "\n";

	return best.costValue;
}

#endif
