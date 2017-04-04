/*
 * (1+1)-ES_PSRS.h
 *
 *  Created on: 04 apr 2017
 *      Author: RedShy
 */

#ifndef SRC__1_1__ES_PSRS_H_
#define SRC__1_1__ES_PSRS_H_

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

bool ES_isValid(ES_MatchingSchema m);

int evolutionStrategy_one_one_psrs(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations, const unsigned maxAttempts,
		const unsigned NThread)
{
	unsigned results[4];

//	for (unsigned i = 0; i < NThread; ++i)
//	{
//		std::thread s(evolutionStrategy_one_one_srs_t, s1, s2, s1l, s2l, sig1, sig2,
//				sig1l, sig2l, p1, m, e, max_generations, maxAttempts, results,
//				i);
//	}

	std::thread t1(evolutionStrategy_one_one_srs_t, s1, s2, s1l, s2l, sig1,
			sig2, sig1l, sig2l, p1, std::ref(m), std::ref(e), max_generations,
			maxAttempts, results, 0);
	std::thread t2(evolutionStrategy_one_one_srs_t, s1, s2, s1l, s2l, sig1,
			sig2, sig1l, sig2l, p1, std::ref(m), std::ref(e), max_generations,
			maxAttempts, results, 1);
	std::thread t3(evolutionStrategy_one_one_srs_t, s1, s2, s1l, s2l, sig1,
			sig2, sig1l, sig2l, p1, std::ref(m), std::ref(e), max_generations,
			maxAttempts, results, 2);
	std::thread t4(evolutionStrategy_one_one_srs_t, s1, s2, s1l, s2l, sig1,
			sig2, sig1l, sig2l, p1, std::ref(m), std::ref(e), max_generations,
			maxAttempts, results, 3);


	t1.join();
	t2.join();
	t3.join();
	t4.join();

	unsigned min = results[0];
	for (unsigned i = 1; i < NThread; i++)
	{
		if (results[i] < min)
		{
			min = results[i];
		}
	}

	return min;
}

void evolutionStrategy_one_one_srs_t(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations, const unsigned maxAttempts,
		unsigned results[], const unsigned index)
{

	const unsigned maxPlateu = 10;
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
	results[index] = best.costValue;
}

bool ES_isValid(ES_MatchingSchema m)
{
	//TODO validate a matching schema
	return true;
}



#endif /* SRC__1_1__ES_PSRS_H_ */
