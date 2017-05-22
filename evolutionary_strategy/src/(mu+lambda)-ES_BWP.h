/*
 * (mu+lambda)-ES_PWP.h
 *
 *  Created on: 04 apr 2017
 *      Author: RedShy
 */

#ifndef mu_lambda_pwp_ES
#define mu_lambda_pwp_ES
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <thread>
#include <queue>
#include "ES_MatchingSchema.h"
#include "EditDistance.h"
#include "MatchingSchema.h"

void evolutionStrategy_WP_t(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations, const unsigned mu,
		const unsigned lambda, unsigned results[], const unsigned index);

int evolutionStrategy_BWP(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations,
		const unsigned mu, const unsigned lambda, const unsigned NThread)
{
	unsigned results[2];

//	for (unsigned i = 0; i < NThread; ++i)
//	{
//		std::thread s(evolutionStrategy_one_one_srs_t, s1, s2, s1l, s2l, sig1, sig2,
//				sig1l, sig2l, p1, m, e, max_generations, maxAttempts, results,
//				i);
//	}

	std::thread t1(evolutionStrategy_WP_t, s1, s2, s1l, s2l, sig1, sig2, sig1l,
			sig2l, p1, std::ref(m), std::ref(e), max_generations,
			15, lambda,
			results, 0);
	std::thread t2(evolutionStrategy_WP_t, s1, s2, s1l, s2l, sig1, sig2, sig1l,
			sig2l, p1, std::ref(m), std::ref(e), max_generations,
			30, lambda,
			results, 1);
//	std::thread t3(evolutionStrategy_WP_t, s1, s2, s1l, s2l, sig1, sig2, sig1l,
//			sig2l, p1, std::ref(m), std::ref(e), max_generations,
//			mu, lambda, results, 2);
//	std::thread t4(evolutionStrategy_WP_t, s1, s2, s1l, s2l, sig1, sig2, sig1l,
//			sig2l, p1, std::ref(m), std::ref(e), max_generations,
//			mu, lambda, results, 3);

	t1.join();
	t2.join();
//	t3.join();
//	t4.join();

	unsigned min = results[0];
	for (unsigned i = 1; i < NThread; i++)
	{
		if (results[i] < min)
		{
			min = results[i];
		}
	}

	std::cout << min;
	return min;
}

void evolutionStrategy_WP_t(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations, const unsigned mu,
		const unsigned lambda, unsigned results[], const unsigned index)
{
	unsigned generation = 0;

	ES_MatchingSchema startingMS(sig1, sig2);

	//Generate mu random individuals
	ES_MatchingSchema parents[mu];
	for (unsigned i = 0; i < mu; ++i)
	{
		startingMS.shuffle();

		//validate matching schema
		if (ES_isValid(startingMS))
		{
			startingMS.costValue = e.edit_distance_matching_schema_enhanced(s1,
					s2, s1l, s2l, startingMS.sigma1, startingMS.sigma2, sig1l,
					sig2l, m);
			parents[i] = startingMS;
		}
		else
		{
			//TODO: not valid, maybe mutate until is valid?
			//repeat iteration
			i--;
		}
		}

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

				//select the worst parent, mu is always very very small like 5 or 10
				unsigned worstParentCostValue = parents[0].costValue;
				unsigned worstParent = 0;
				for (unsigned i = 1; i < mu; i++)
					{
					if (parents[i].costValue > worstParentCostValue)
					{
						worstParentCostValue = parents[i].costValue;
						worstParent = i;
					}
					}

				int newDistance =
						e.edit_distance_matching_schema_enhanced_with_diagonal(
								s1, s2, s1l, s2l, child.sigma1, child.sigma2,
								sig1l, sig2l, m, worstParentCostValue);

				if (newDistance != -1)
				{
					//The child is better than the worst parent, so he become a new parent
					child.costValue = newDistance;
					parents[worstParent] = child;
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

	//TODO return best of all
	unsigned bestValue = parents[0].costValue;
	unsigned bestParent = 0;
	for (unsigned i = 1; i < mu; i++)
	{
		if (parents[i].costValue < bestValue)
			{
			bestValue = parents[i].costValue;
			bestParent = i;
			}
		}

	results[index] = bestValue;
}

#endif
