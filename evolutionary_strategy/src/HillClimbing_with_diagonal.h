/*
 * HillClimbing_with_diagonal.h
 *
 *  Created on: 05 apr 2017
 *      Author: RedShy
 */

#ifndef SRC_HILLCLIMBING_WITH_DIAGONAL_H_
#define SRC_HILLCLIMBING_WITH_DIAGONAL_H_

#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <ctime>
#include <vector>
#include <bitset>
#include "Alignment.h"
#include "EditDistance.h"
#include "FixedED.h"
#include "MatchingSchema.h"
#include "ES_MatchingSchema.h"
#include "Utility.h"

#define CLOCKS_PER_MS (CLOCKS_PER_SEC / 1000)

int hill_climbing_d(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,
		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l, const size_t& p1,
		matching_schema<bool>& m, edit_distance& e)
{
	clock_t start = clock();
	long double msElapsed = 0;
	//std::cout << "enter the void (1)" << endl;

	unsigned d = e.edit_distance_matching_schema(s1, s2, s1l, s2l, m);
	unsigned minDist = d;
	unsigned minMinDist = minDist;

	//std::cout << "enter the void (2)" << endl;

	// for the permutations
	unsigned* sigma1_o = new unsigned[sig1l];
	std::iota(sigma1_o, sigma1_o + sig1l, 0);
	unsigned* sigma2_o = new unsigned[sig2l];
	std::iota(sigma2_o, sigma2_o + sig2l, 0);
	unsigned* sigma1_t = new unsigned[sig1l];
	std::iota(sigma1_t, sigma1_t + sig1l, 0);
	unsigned* sigma2_t = new unsigned[sig2l];
	std::iota(sigma2_t, sigma2_t + sig2l, 0);

	//std::cout << "enter the void (3)" << endl;

	// for fixpoints
	unsigned* sigma1_min = new unsigned[sig1l];
	std::iota(sigma1_min, sigma1_min + sig1l, 0);
	unsigned* sigma2_min = new unsigned[sig2l];
	std::iota(sigma2_min, sigma2_min + sig2l, 0);
	unsigned* sigma1_min_min = new unsigned[sig1l];
	std::iota(sigma1_min_min, sigma1_min_min + sig1l, 0);
	unsigned* sigma2_min_min = new unsigned[sig2l];
	std::iota(sigma2_min_min, sigma2_min_min + sig2l, 0);

	//std::cout << "enter the void (4)" << endl;

	size_t attempts = 1, shuffle_tries = 2;
	unsigned tries = 0, k_shuffle = 0;

	//std::cout << "start: " << d << endl;

	clock_t timeElapsed = clock() - start;
	msElapsed = timeElapsed / CLOCKS_PER_MS;
	std::cout << msElapsed << " " << minMinDist << "\n";

	bool improved = true;
	while (improved)
	{
		improved = false;

		for (size_t ip = 0; ip < sig1l; ip++)
		{
			for (size_t jp = ip; jp < sig1l; jp++)
			{

				std::copy(sigma1_t, sigma1_t + sig1l, sigma1_o); // reset state
				std::swap(sigma1_o[ip], sigma1_o[jp]);					// swap

				if (isValid(sigma1_o, sig1l, p1))
				{

					for (size_t ipp = 0; ipp < sig2l; ipp++)
					{
						for (size_t jpp = ipp; jpp < sig2l; jpp++)
						{

							std::copy(sigma2_t, sigma2_t + sig2l, sigma2_o);// reset state
							std::swap(sigma2_o[ipp], sigma2_o[jpp]);	// swap

							int newDistance =
									e.edit_distance_matching_schema_enhanced_with_diagonal(
											s1, s2, s1l, s2l, sigma1_o,
											sigma2_o, sig1l, sig2l, m, minDist);

							if (newDistance != -1)
							{
								minDist = newDistance;

								improved = true;
								std::copy(sigma1_o, sigma1_o + sig1l,
										sigma1_min);
								std::copy(sigma2_o, sigma2_o + sig2l,
										sigma2_min);

							}
						}
						std::copy(sigma2_t, sigma2_t + sig2l, sigma2_o);
					}

				}
			}

		}

		if (improved)
		{
			// copy sigmaMin to sigmaOrig
			std::copy(sigma1_min, sigma1_min + sig1l, sigma1_o);
			std::copy(sigma2_min, sigma2_min + sig2l, sigma2_o);
			// copy sigmaOrig to sigmaTmp
			std::copy(sigma1_o, sigma1_o + sig1l, sigma1_t);
			std::copy(sigma2_o, sigma2_o + sig2l, sigma2_t);
		}
		else
		{
			if (minDist < minMinDist)
			{
				minMinDist = minDist;

				// copy sigmaMin to sigmaMinMin
				std::copy(sigma1_min, sigma1_min + sig1l, sigma1_min_min);
				std::copy(sigma2_min, sigma2_min + sig2l, sigma2_min_min);

				improved = true;
				tries = 0;
			}

			if (tries < attempts)
			{
				improved = true;
				tries++;

				// random swap
				// for sigma1, we try _SHUFFLE_TRIES times, then if is still not valid, we retry with the original one
				for (k_shuffle = 0;
						k_shuffle < shuffle_tries
								&& !isValid(sigma1_t, sig1l, p1); ++k_shuffle)
					shuffle(sigma1_t, sig1l);
				if (k_shuffle == shuffle_tries) std::copy(sigma1_o,
						sigma1_o + sig1l, sigma1_t);
				// no constraints on the shuffle for sigma2
				shuffle(sigma2_t, sig2l);

				std::copy(sigma1_t, sigma1_t + sig1l, sigma1_o);
				std::copy(sigma2_t, sigma2_t + sig2l, sigma2_o);

				minDist = e.edit_distance_matching_schema_enhanced(s1, s2, s1l,
						s2l, sigma1_o, sigma2_o, sig1l, sig2l, m);
			}
		}
		clock_t timeElapsed = clock() - start;
		msElapsed = timeElapsed / CLOCKS_PER_MS;
		std::cout << msElapsed << " " << minMinDist << "\n";

	}

	delete[] sigma1_o;
	delete[] sigma1_t;
	delete[] sigma2_o;
	delete[] sigma2_t;
	delete[] sigma1_min;
	delete[] sigma1_min_min;
	delete[] sigma2_min;
	delete[] sigma2_min_min;

	return minMinDist;
}


#endif /* SRC_HILLCLIMBING_WITH_DIAGONAL_H_ */
