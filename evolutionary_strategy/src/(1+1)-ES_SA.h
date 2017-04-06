/*
 * (1+1)-ES_SA.h
 *
 *  Created on: 05 apr 2017
 *      Author: RedShy
 */

#ifndef SRC__1_1__ES_SA_H_
#define SRC__1_1__ES_SA_H_


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

int evolutionStrategy_one_one_sa(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations)
{
	unsigned generation = 0;

	const unsigned G = 10;
	const double s = (1.00 / 5.00);
	unsigned G_s = 0;
	unsigned mutationStrength = 1;

	unsigned plateu = 0;
	const unsigned maxPlateu = 40;

	ES_MatchingSchema parent(sig1, sig2);
	//Random start
	parent.shuffle();

	parent.costValue = e.edit_distance_matching_schema_enhanced(s1, s2, s1l,
			s2l, parent.sigma1, parent.sigma2, sig1l, sig2l, m);
	while (generation <= max_generations)
	{
		//Produce child
		ES_MatchingSchema child = parent;

		//mutate child
		child.mutate(mutationStrength);

		//validate child
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
				G_s++;
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

		if (generation % G == 0)
		{

			const float r = G_s / G;
			std::cout << "DENTRO! G_s=" << G_s << " G=" << G << " r=" << r
					<< " s=" << s << "\n";
			if (r > s)
			{
				std::cout << "INCREMENTO!\n";
				++mutationStrength;
			}
			else if (r < s)
			{
				std::cout << "VORREI DECREMENTARE\n";
				if (mutationStrength != 1)
				{
					std::cout << "DECREMENTO!\n";
					--mutationStrength;
				}
			}
			G_s = 0;
		}
	}

	//TODO return best of all
	return parent.costValue;
}




#endif /* SRC__1_1__ES_SA_H_ */
