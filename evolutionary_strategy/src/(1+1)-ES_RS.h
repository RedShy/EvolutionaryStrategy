/*
 * (1+1)-ES.h
 *
 *  Created on: 03 apr 2017
 *      Author: RedShy
 */

#ifndef ES_ONE_ONE_RS
#define ES_ONE_ONE_RS

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

int evolutionStrategy_one_one_rs(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations)
{
	clock_t start = clock();
	long double msElapsed = 0;

	unsigned generation = 0;
	unsigned plateu = 0;
	unsigned lastBest = std::numeric_limits<unsigned int>::max();
	const unsigned threshold = 0;
	const unsigned lapCheckPoint = 500;
	const unsigned maxPlateu = 50 * p1;



	ES_MatchingSchema parent(sig1, sig2);

	//Random start
	parent.shuffle();
	parent.costValue = e.edit_distance_matching_schema_enhanced(s1, s2, s1l,
			s2l, parent.sigma1, parent.sigma2, sig1l, sig2l, m);

	ES_MatchingSchema best = parent;

	clock_t timeElapsed1 = clock() - start;
	msElapsed = timeElapsed1 / CLOCKS_PER_MS;
	std::cout << msElapsed << " " << parent.costValue << "\n";

	while (generation <= max_generations)
	{
		//Produce child
		ES_MatchingSchema child = parent;

		//mutate child
		child.mutate();

		if (ES_isValid(child))
		{
			int newDistance =
					e.edit_distance_matching_schema_enhanced_with_diagonal(s1,
							s2, s1l, s2l, child.sigma1, child.sigma2, sig1l,
							sig2l, m, parent.costValue);

			if (newDistance != -1)
			{
				//The child is better than its father, so he become new parent
				parent = child;
				parent.costValue = newDistance;

				plateu = 0;

				//TODO maybe we can do better
				if (parent.costValue < best.costValue)
				{
					best = parent;

					clock_t timeElapsed = clock() - start;
					msElapsed = timeElapsed / CLOCKS_PER_MS;
					std::cout << msElapsed << " " << parent.costValue << "\n";

				}
			}
			else
			{
				plateu++;
				if (plateu == maxPlateu)
				{
					plateu = 0;

					parent.shuffle();
					parent.costValue = e.edit_distance_matching_schema_enhanced(
							s1, s2, s1l, s2l, parent.sigma1, parent.sigma2,
							sig1l, sig2l, m);
					if (parent.costValue < best.costValue)
					{
						best = parent;

						clock_t timeElapsed = clock() - start;
						msElapsed = timeElapsed / CLOCKS_PER_MS;
						std::cout << msElapsed << " " << parent.costValue
								<< "\n";

					}
				}
			}
		}
		else
		{
			continue;
		}

		generation++;

		//TODO experimental return on investiment
		if (generation % lapCheckPoint == 0)
		{
			//few rendiment in the last checkpoint, so stop here
			if (lastBest - best.costValue < threshold)
			{
				break;
			}

			//new value for checkpoint
			lastBest = best.costValue;
		}

	}

	//TODO return best of all
	clock_t timeElapsed = clock() - start;
	msElapsed = timeElapsed / CLOCKS_PER_MS;
	std::cout << msElapsed << " " << best.costValue << "\n";

	return best.costValue;
}




#endif /* SRC__1_1__ES_H_ */
