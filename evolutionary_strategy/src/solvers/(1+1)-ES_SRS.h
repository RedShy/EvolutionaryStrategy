
/*
 * (1+1)-ES_RA.h
 *
 *  Created on: 04 apr 2017
 *      Author: RedShy
 */
#ifndef ES_ONE_ONE_SRS
#define ES_ONE_ONE_SRS

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

int evolutionStrategy_one_one_srs(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations, const unsigned maxAttempts)
{

	const unsigned maxPlateu = 5 * p1;
	unsigned attempts = 0;
	ES_MatchingSchema parent(sig1, sig2);
	//Random start
	parent.shuffle();
	parent.costValue = e.edit_distance_matching_schema_enhanced(s1, s2, s1l,
			s2l, parent.sigma1, parent.sigma2, sig1l, sig2l, m);
	ES_MatchingSchema best = parent;

	while (attempts < maxAttempts)
	{
		unsigned generation = 0;
		unsigned plateu = 0;

		while (generation <= max_generations)
		{
			//Produce child
			ES_MatchingSchema child = parent;

			//mutate child
			child.mutate();

			//validate child
			if (ES_isValid(child))
			{
				int newDistance =
						e.edit_distance_matching_schema_enhanced_with_diagonal(
								s1, s2, s1l, s2l, child.sigma1, child.sigma2,
								sig1l, sig2l, m, parent.costValue);

				if (newDistance != -1)
				{
					//The child is better than its father, so he become new parent
					parent = child;
					parent.costValue = newDistance;

					plateu = 0;
				}
				else
				{
					plateu++;
					if (plateu == maxPlateu)
					{
						break;
					}
				}
				//else the child is worse than its father so he is discarded
			}
			else
			{
				//TODO: not valid, maybe mutate until is valid?
				//repeat iteration
				continue;
			}

			generation++;
		}

		//check if the last attempt has improved the solution
		if (parent.costValue < best.costValue)
		{
			best = parent;
		}

		//Random restart
		parent.shuffle();
		parent.costValue = e.edit_distance_matching_schema_enhanced(s1, s2, s1l,
				s2l, parent.sigma1, parent.sigma2, sig1l, sig2l, m);
		attempts++;
	}
	return best.costValue;
}

#endif /* SRC__1_1__ES_H_ */

