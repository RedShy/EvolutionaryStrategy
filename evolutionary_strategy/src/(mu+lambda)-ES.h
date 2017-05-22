#ifndef mu_lambda_ES
#define mu_lambda_ES

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


std::vector<ES_MatchingSchema> selectBestIndividuals(const unsigned mu,
		std::vector<ES_MatchingSchema>& individuals, const bool plusSelection);

int evolutionStrategy(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations, const unsigned mu,
		const unsigned lambda, const bool plusSelection)
{
	unsigned generation = 0;

	ES_MatchingSchema startingMS(sig1, sig2);

	//Generate mu random individuals
	std::vector<ES_MatchingSchema> parents;
	for (unsigned i = 0; i < mu; ++i)
	{
		startingMS.shuffle();

		//validate matching schema
		if (ES_isValid(startingMS))
		{
			startingMS.costValue =
					e.edit_distance_matching_schema_enhanced(s1,
							s2, s1l, s2l, startingMS.sigma1, startingMS.sigma2,
							sig1l, sig2l, m);
			parents.push_back(startingMS);
		}
		else
		{
			//TODO: not valid, maybe mutate until is valid?
			//repeat iteration
			i--;
		}
	}

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

				child.costValue =
						e.edit_distance_matching_schema_enhanced(
								s1, s2, s1l, s2l, child.sigma1, child.sigma2,
								sig1l, sig2l, m);

				parents.push_back(child);
			}
			else
			{
				//TODO: not valid, maybe mutate until is valid?
				//repeat iteration
				i--;
			}
		}


		parents = selectBestIndividuals(mu, parents, plusSelection);


		generation++;
	}

	//TODO return best of all
	make_heap(parents.begin(), parents.end());
	return parents.front().costValue;
}

std::vector<ES_MatchingSchema> selectBestIndividuals(const unsigned mu,
		std::vector<ES_MatchingSchema>& individuals, const bool plusSelection)
{
	//TODO using a heap
	std::vector<ES_MatchingSchema> bestIndividuals;
	if (plusSelection)
	{
		make_heap(individuals.begin(), individuals.end());

		for (unsigned i = 0; i < mu; i++)
		{
			pop_heap(individuals.begin(), individuals.end());
			bestIndividuals.push_back(individuals.back());
			individuals.pop_back();
		}
	}
	else
	{
		//TODO Comma selection
	}

	return bestIndividuals;
}
#endif
