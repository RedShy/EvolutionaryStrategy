/*
 * (mu+lambda)-ES-comma.h
 *
 *  Created on: 20 giu 2017
 *      Author: RedShy
 */

#ifndef SRC__MU_LAMBDA__ES_COMMA_H_
#define SRC__MU_LAMBDA__ES_COMMA_H_

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

int evolutionStrategy_comma(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, const size_t& p2, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations, const unsigned mu,
		const unsigned lambda)
{
	if(lambda < mu)
	{
		std::cout<<"lambda has to be greater than mu\n";
		exit(1);
	}

	clock_t start = clock();
	long double msElapsed = 0;

	//Initialize stuff for the mutator swap2-E
	const unsigned * const blocksig1 = initializeBlocksSwap2E(sig1, p1);
	const unsigned * const blocksig2 = initializeBlocksSwap2E(sig2, p2);

	ES_MatchingSchema startingMS(sig1, sig2);

	ES_MatchingSchema best;
	best.costValue = std::numeric_limits<unsigned int>::max();

	//Generate mu random individuals
	ES_MatchingSchema* parents = new ES_MatchingSchema[mu];
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
	const unsigned remainingChildren = lambda-mu;


	ES_MatchingSchema* children = new ES_MatchingSchema[mu];
	unsigned generation = 0;
	while (generation <= max_generations)
	{
		for (unsigned i = 0; i < mu; i++)
		{
			//Choose random parent
			const unsigned p = rand() % mu;

			//Produce child, in the case parents=1 (like this) just clone
			children[i] = parents[p];

			//mutate child
			children[i].swap2_enhanced(blocksig1, blocksig2);

			children[i].costValue =
					e.edit_distance_matching_schema_enhanced(s1,
							s2, s1l, s2l, children[i].sigma1, children[i].sigma2, sig1l,
							sig2l, m);

			if (children[i].costValue < best.costValue)
			{
				best = children[i];

				clock_t timeElapsed = clock() - start;
				msElapsed = timeElapsed / CLOCKS_PER_MS;
				std::cout << msElapsed << " " << best.costValue << "\n";
			}
		}

		//Generate the remaining children. Only mutation, no recombination
		for (unsigned i = 0; i < remainingChildren; i++)
		{
			//select the worst child
			unsigned worstChild=0;
			unsigned worstChildCostValue=children[0].costValue;
			for(unsigned i=1; i<mu; ++i)
			{
				if(children[i].costValue > worstChildCostValue)
				{
					worstChildCostValue=children[i].costValue;
					worstChild=i;
				}
			}

			//Choose random parent
			const unsigned p = rand() % mu;

			//Produce child, in the case parents=1 (like this) just clone
			ES_MatchingSchema child = parents[p];

			//mutate child
			child.swap2_enhanced(blocksig1, blocksig2);

			const int newDistance =
					e.edit_distance_matching_schema_enhanced_with_diagonal(s1,
							s2, s1l, s2l, child.sigma1, child.sigma2, sig1l,
							sig2l, m, worstChildCostValue);

			if (newDistance != -1)
			{
				//The child is better than the worst child,
				child.costValue = newDistance;

				//so he replace the worst child
				children[worstChild] = child;

				if (child.costValue < best.costValue)
				{
					best = child;

					clock_t timeElapsed = clock() - start;
					msElapsed = timeElapsed / CLOCKS_PER_MS;
					std::cout << msElapsed << " " << best.costValue << "\n";
				}
			}
		}

		ES_MatchingSchema* tmp=children;
		children=parents;
		parents=tmp;

		generation++;
	}

	delete[] blocksig1;
	delete[] blocksig2;

	delete[] parents;
	delete[] children;

	clock_t timeElapsed = clock() - start;
	msElapsed = timeElapsed / CLOCKS_PER_MS;
	std::cout << msElapsed << " " << best.costValue << "\n";

//	m.print_matching_schema(best.sigma1,best.sigma2);

//	std::cout<<best.costValue;
//	std::cout<<"CHIAMATE A EDIT DISTANCE= "<<edit_distance::tentativi<<"\n";
	return best.costValue;
}



#endif /* SRC__MU_LAMBDA__ES_COMMA_H_ */
