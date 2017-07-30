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
#include <map>

#define CLOCKS_PER_MS (CLOCKS_PER_SEC / 1000)

int bruteforce(const std::vector<unsigned>& s1, const std::vector<unsigned>& s2,
		const size_t& s1l, const size_t& s2l, const std::vector<unsigned>& sig1,
		const std::vector<unsigned>& sig2, const size_t& sig1l,
		const size_t& sig2l, const matching_schema<bool>& m, edit_distance& e)
{

	clock_t start = clock();
	long double msElapsed = 0;

	unsigned distance = e.edit_distance_matching_schema(s1, s2, s1l, s2l, m);
	unsigned current = distance;

	unsigned* perm1 = new unsigned[sig1l];
	for (unsigned i = 0; i < sig1l; ++i)
		perm1[i] = i;
	unsigned* perm2 = new unsigned[sig2l];
	for (unsigned i = 0; i < sig2l; ++i)
		perm2[i] = i;

	FixedED<unsigned> fixed_ed(s1l + 1, s2l + 1);

	unsigned int computed=1;
	unsigned int totalPermutations=1;
	for(int i=1; i<=sig1l; i++)
	{
		totalPermutations=i*totalPermutations;

	}

	std::map<int,int> editDistances=std::map<int,int>();
	do
	{
//		do
//		{
			current = e.edit_distance_matching_schema_enhanced(s1, s2,
					s1l, s2l, perm1, perm2, sig1l, sig2l, m);

//			if(current!=-1)
//			{
//				distance = current;
//			}


			if(current < distance)
			{
				distance=current;
			}

			editDistances[current]++;

			computed++;

			if(computed % 1000 == 0)
			{
				clock_t timeElapsed = clock() - start;
				msElapsed = timeElapsed / CLOCKS_PER_MS;

				std::cout << " Computing permutation number " << computed << " on "<<totalPermutations <<" total. Edit Distance: "<<distance<<" Time: "<<msElapsed<< "\n";
			}

		} while (std::next_permutation(perm2, perm2 + sig2l));

//	} while (std::next_permutation(perm1, perm1 + sig1l));

	std::cout<<"RESULTS:"<<std::endl;

	for(std::map<int,int>::iterator i=editDistances.begin(); i!=editDistances.end(); ++i)
	{
		double percent= (i->second/(double) totalPermutations) * 100;
		std::cout<<"Edit distance: "<<i->first <<" Encountered: "<<i->second<<" Percentage: "<<percent<<" % of total"<<std::endl;
	}

//	clock_t timeElapsed = clock() - start;
//	msElapsed = timeElapsed / CLOCKS_PER_MS;
//	std::cout << msElapsed << " " << distance << "\n";

	return distance;
}



#endif /* SRC_BRUTEFORCE_H_ */
