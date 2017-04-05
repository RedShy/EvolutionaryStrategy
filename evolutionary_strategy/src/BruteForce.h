/*
 * BruteForce.h
 *
 *  Created on: 05 apr 2017
 *      Author: HantolR
 */

#ifndef SRC_BRUTEFORCE_H_
#define SRC_BRUTEFORCE_H_
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <ctime>
#include <vector>
#include <bitset>

int bruteforce(const std::vector<unsigned>& s1, const std::vector<unsigned>& s2,
		const size_t& s1l, const size_t& s2l, const std::vector<unsigned>& sig1,
		const std::vector<unsigned>& sig2, const size_t& sig1l,
		const size_t& sig2l, const matching_schema<bool>& m, edit_distance& e)
{

	unsigned distance = e.edit_distance_matching_schema(s1, s2, s1l, s2l, m);
	unsigned current = distance;

	unsigned* perm1 = new unsigned[sig1l];
	for (unsigned i = 0; i < sig1l; ++i)
		perm1[i] = i;
	unsigned* perm2 = new unsigned[sig2l];
	for (unsigned i = 0; i < sig2l; ++i)
		perm2[i] = i;

	FixedED<unsigned> fixed_ed(s1l + 1, s2l + 1);

	do
	{
		do
		{
			current = fixed_ed.edit_distance_matching_schema_enhanced(s1, s2,
					s1l, s2l, perm1, perm2, sig1l, sig2l, m);

			if (current < distance) distance = current;

		} while (std::next_permutation(perm2, perm2 + sig2l));
	} while (std::next_permutation(perm1, perm1 + sig1l));

	return distance;
}



#endif /* SRC_BRUTEFORCE_H_ */
