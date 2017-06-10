#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <ctime>
#include <vector>
#include <bitset>
#include <thread>
#include <time.h>
#include "Alignment.h"
#include "EditDistance.h"
#include "FixedED.h"
#include "MatchingSchema.h"
#include "ES_MatchingSchema.h"
#include "Utility.h"
/* Solvers */
#include "HillClimbing.h"
#include "BruteForce.h"
#include "(mu+lambda)-ES_WP_RS.h"
#include "(1+1)-ES.h"
#include "(1+1)-ES_SRS.h"
#include "(mu+lambda)-ES_AF.h"
#include "(mu+lambda)-ES.h"
#include "(mu+lambda)-ES-shuffle.h"
#include "(1+1)-ES_RS.h"
#include "(mu+1)-ES_WP.h"
#include "random_search.h"
#include "swap2-2.h"
#include "swap2-3.h"
#include "swap2-4.h"

/* Definitions */
#define endl '\n'
#define CLOCKS_PER_MS (CLOCKS_PER_SEC / 1000)

/* Consts */
const unsigned short _ASCII_LEN = 255 - 0;
const bool _DEBUG = false;
const std::string _HC_ARG("hc");
const std::string _BRUTEFORCE_ARG("ex");
const std::string _ES_WP_ARG("es-wp");
const std::string _ES_ARG("es");
const std::string _ES_ONE_ONE_ARG("es_one_one");
const std::string _ES_ONE_ONE_SRS_ARG("es_one_one_srs");
const std::string _ES_AF_ARG("es_af");
const std::string _ES_ONE_ONE_SA_ARG("es_one_one_sa");
const std::string _ES_ONE_ONE_RS_ARG("es_one_one_rs");
const std::string _SPECIFIC_PMS("specific-permutations");
const std::string _SPECIFIC_MMS("specific-matrix");

/* Utility functions */
void extract_sigma(std::string&, std::string&);
void define_mapping(std::string&, std::map<char, int>&);
void extract_sigma_fibre(std::vector<std::string>&, std::vector<std::string>&);
void define_mapping_fibre(std::vector<std::string>&,
		std::map<std::string, int>&);
void print_alignment(const Alignment<int>&, const std::string&,
		const std::string&, const matching_schema<bool>&, int&, bool);



int main(int argc, char *argv[])
{
//	std::ios_base::sync_with_stdio(false);
	srand(time(0));

	// arguments: [hc|ex] p1 p2 [specific-permutations|specific-matrix]
	// heuristic
	std::string heuristic(argv[1]);
	// p1 and p2
	size_t p1 = argc > 3 ? fast_atoi(argv[2]) : 1;
	size_t p2 = argc > 3 ? fast_atoi(argv[3]) : 1;
	// specific matching schema
	bool specific_perm = false;
	bool specific_matrix = false;
	if (argc > 4)
	{
		std::string specific_temp(argv[4]);
		// if we have in input a specific ms with permutations
		if (specific_temp == _SPECIFIC_PMS) specific_perm = true;
		// if we have in input a specific ms with the whole matrix
		if (specific_temp == _SPECIFIC_MMS) specific_matrix = true;
	}

	// TODO: implement constraints managing
	bool has_constraints = false;
	bool default_constraints_mode = false;

	/* (1) read strings, extract sigmas and constraints*/
	std::string s1, s2;
	std::string sigma1(""), sigma2("");

	read_stdin(s1, s2);
	extract_sigma(s1, sigma1);
	extract_sigma(s2, sigma2);

	size_t sigma1l = sigma1.size();
	size_t sigma2l = sigma2.size();
	size_t s1l = s1.size();
	size_t s2l = s2.size();

	std::vector<p_constr> constraints;
	if (has_constraints) read_constraints(constraints);

	/* define the mapping from char -> int */
	std::map<char, int> map1;
	std::map<char, int> map2;

	define_mapping(sigma1, map1);
	define_mapping(sigma2, map2);

	/* integer representations of strings and sigmas */
	std::vector<unsigned> s1i(s1l);
	std::vector<unsigned> s2i(s1l);
	std::vector<unsigned> sigma1i(sigma1l);
	std::vector<unsigned> sigma2i(sigma2l);

	std::iota(sigma1i.begin(), sigma1i.end(), 0);		// encoding for sigma1
	std::iota(sigma2i.begin(), sigma2i.end(), 0);		// encoding for sigma2
	for (size_t i = 0; i < s1l; ++i)
		s1i[i] = map1[s1[i]];		// encoding for s1
	for (size_t i = 0; i < s2l; ++i)
		s2i[i] = map2[s2[i]];		// encoding for s2

	if (_DEBUG)
	{
		std::cout << "p1: " << p1 << ", p2: " << p2 << endl;
		std::cout << "s len: " << s1.size() << ", " << s2.size() << endl;
		std::cout << "strings: " << s1 << ", " << s2 << endl;
		std::cout << "sigma len: " << sigma1.size() << ", " << sigma2.size()
				<< endl;
		std::cout << "sigmas: " << sigma1 << ", " << sigma2 << endl;
		std::cout << "int rep (s1): ";
		for (size_t i = 0; i < s1l; ++i)
			std::cout << s1i[i];
		std::cout << endl;
		std::cout << "int rep (s2): ";
		for (size_t i = 0; i < s2l; ++i)
			std::cout << s2i[i];
		std::cout << endl;
		std::cout << "int rep (sigma1): ";
		for (size_t i = 0; i < sigma1l; ++i)
			std::cout << sigma1i[i];
		std::cout << endl;
		std::cout << "int rep (sigma2): ";
		for (size_t i = 0; i < sigma2l; ++i)
			std::cout << sigma2i[i];
		std::cout << endl;
	}

	/* identity (classical) matching schema */
	matching_schema<bool> ms(sigma1l, sigma2l, p1, p2, true,
			default_constraints_mode);
	ms.set_general(sigma1, sigma2, false);



	// TODO: CONSTRAINTS
	//if (has_constraints)
	//	ms.set_constraints(map1, map2, constraints, !default_constraints_mode);

	if (_DEBUG)
	{
		if (has_constraints)
		{
			std::cout << "Constraints: ";
			for (std::vector<p_constr>::size_type i = 0; i < constraints.size();
					++i)
				std::cout << constraints[i].first << ", "
						<< constraints[i].second << endl;
		}
		ms.print_matching_schema(sigma1, sigma2);
	}

	/* create an edit distance object */
	edit_distance e;

	struct timespec start1, finish1;
	double elapsed;

	/* here we call the solver needed */
	int distance = -1;
	double msElapsed = 0;
	// Common execution of HC or EX
	if (!specific_perm && !specific_matrix)
		{
		clock_gettime(CLOCK_MONOTONIC, &start1);
		clock_t start = clock();
		if (heuristic == _BRUTEFORCE_ARG)
		{
			distance = bruteforce(s1i, s2i, s1l, s2l, sigma1i, sigma2i, sigma1l,
					sigma2l, ms, e);
		}
		else if (heuristic == _HC_ARG)
		{
			distance = hill_climbing(s1i, s2i, s1l, s2l, sigma1i, sigma2i,
					sigma1l, sigma2l, p1, ms, e);
		}
		else if (heuristic == _ES_WP_ARG)
		{
			distance = evolutionStrategy_WP(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 14400, 10);
		}
		else if (heuristic == _ES_ARG)
		{
			distance = evolutionStrategy(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 129, 30, 111);
		}
		else if (heuristic == "es-shuffle")
		{
			distance = evolutionStrategy_shuffle(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 200, 40, 200);
		}
		else if (heuristic == "es-400-10-18")
		{
			distance = evolutionStrategy(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 400, 10, 18);
		}
		else if (heuristic == "es-85-10-85")
		{
			distance = evolutionStrategy(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 85, 10, 85);
		}
		else if (heuristic == "es-wp-rs")
		{
			distance = evolutionStrategy_WP_RS(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 57600, 20, 1);
		}
		else if (heuristic == "es-wp-S")
		{
			unsigned sum=0;
			for(int i=0; i<1; i++)
			{
				sum += evolutionStrategy_WP(s1i, s2i, s1l, s2l, sigma1i,
						sigma2i, sigma1l, sigma2l, p1, ms, e, 5000, 10);
			}
			distance=sum/10;
			std::cout<<distance;
		}
		else if (heuristic == "es-wp-M")
		{
			unsigned sum=0;
			for(int i=0; i<1; i++)
			{
				sum += evolutionStrategy_WP(s1i, s2i, s1l, s2l, sigma1i,
						sigma2i, sigma1l, sigma2l, p1, ms, e, 2500, 10);
			}
			distance=sum/10;
			std::cout<<distance;
		}
		else if (heuristic == "es-wp-L")
		{
			unsigned sum=0;
			for(int i=0; i<1; i++)
			{
				sum += evolutionStrategy_WP(s1i, s2i, s1l, s2l, sigma1i,
						sigma2i, sigma1l, sigma2l, p1, ms, e, 1000, 10);
			}
			distance=sum/10;
			std::cout<<distance;
		}
		else if (heuristic == "swap2-2")
		{
			distance = evolutionStrategy_WP_swap2_2(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 50, 5, 72);
			std::cout<<distance;
		}
		else if (heuristic == "swap2-3")
		{
			distance = evolutionStrategy_WP_swap2_3(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 50, 5, 72);
			std::cout<<distance;
		}
		else if (heuristic == "swap2-4")
		{
			distance = evolutionStrategy_WP_swap2_4(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 50, 5, 72);

		}
		else if (heuristic == "es-wp-7200-10-1")
		{
			distance = evolutionStrategy_WP(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 7200, 10);
		}
		else if (heuristic == "es-wp-30000")
		{
			distance = evolutionStrategy_WP_RS(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 30000, 10, 1);
		}
		else if (heuristic == "es-one-one")
		{
			distance = evolutionStrategy_one_one(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 7200, 0);
		}
		else if (heuristic == "es-one-one-srs")
		{
			distance = evolutionStrategy_one_one_srs(s1i, s2i, s1l, s2l,
					sigma1i, sigma2i, sigma1l, sigma2l, p1, ms, e, 1440, 10);
		}
		else if (heuristic == _ES_ONE_ONE_ARG)
		{
			distance = evolutionStrategy_one_one(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 7200, 0);
		}
		else if (heuristic == _ES_ONE_ONE_SRS_ARG)
		{
			distance = evolutionStrategy_one_one_srs(s1i, s2i, s1l, s2l,
					sigma1i, sigma2i, sigma1l, sigma2l, p1, ms, e, 1440, 5);
		}
		else if (heuristic == _ES_AF_ARG)
		{
			distance = evolutionStrategy_AF(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 20, 5, 18);
		}
		else if (heuristic == _ES_ONE_ONE_RS_ARG)
		{
			distance = evolutionStrategy_one_one_rs(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 5000);
		}
		else if (heuristic == "random")
		{
			distance = random_search(s1i, s2i, s1l, s2l, sigma1i,
					sigma2i, sigma1l, sigma2l, p1, ms, e, 2000000);
		}
		clock_t timeElapsed = clock() - start;
		clock_gettime(CLOCK_MONOTONIC, &finish1);
		msElapsed = timeElapsed / CLOCKS_PER_MS;

		elapsed = (finish1.tv_sec - start1.tv_sec);
		elapsed += (finish1.tv_nsec - start1.tv_nsec) / 1000000000.0;

	}
	// For a specific matching schema
	else
	{
		// input matching schema
		unsigned* sigma1_spe = new unsigned[sigma1l];
		unsigned* sigma2_spe = new unsigned[sigma2l];

		// if we have the permutations
		if (specific_perm)
		{
			// read string
			std::string read;

			// skip the obtained edit distance
			getline(std::cin, read);

			// read the specified matching schema
			read_specific_matchingschema(read, sigma1_spe);
			read_specific_matchingschema(read, sigma2_spe);

			// check if there is a matching schema or not
			bool there_is_a_ms = false;
			for (size_t i = 0; i < sigma1l; ++i)
				if (sigma1_spe[i])
				{
					there_is_a_ms = true;
					break;
				}

			if (there_is_a_ms)
			{
				distance = e.edit_distance_matching_schema_enhanced(s1i, s2i,
						s1l, s2l, sigma1_spe, sigma2_spe, sigma1l, sigma2l, ms);
			}
			else
				distance = e.edit_distance_matching_schema(s1i, s2i, s1l, s2l,
						ms);

			// if we have the matrix
		}
		else if (specific_matrix)
		{
			// populate linearly
			for (size_t i = 0; i < sigma1l; ++i)
				sigma1_spe[i] = i;
			for (size_t i = 0; i < sigma2l; ++i)
				sigma2_spe[i] = i;

			// read the matrix
			read_specific_matrix(ms);

			// and now for something completely different...
			distance = e.edit_distance_matching_schema_enhanced(s1i, s2i, s1l,
					s2l, sigma1_spe, sigma2_spe, sigma1l, sigma2l, ms);
		}
	}

//	std::cout << distance;
//	if (heuristic == _ES_ONE_ONE_BSRS_ARG || heuristic == _ES_BWP_ARG)
//	{
//		std::cout << ' ' << (int) (elapsed * 1000) << endl;
//	}
//	else
//	{
//	std::cout << ' ' << msElapsed << endl;
//	}

	return 0;
}

/**
 * extract_sigmas works with a bitset of length _ASCII_LEN = 255 - _ASCII_START
 */
void extract_sigma(std::string& s, std::string& e)
{
	std::bitset<_ASCII_LEN> symbols;

	for (size_t i = 0; i < s.size(); ++i)
		symbols[(int) s[i]] = (symbols[(int) s[i]] || 1);

	for (size_t i = 0; i < _ASCII_LEN; ++i)
		if (symbols[i]) e += (char) i;
}

void extract_sigma_fibre(std::vector<std::string>& s,
		std::vector<std::string>& e)
{
	for (int i = 0; i < s.size(); i++)
		e.push_back(s[i]);

	sort(e.begin(), e.end());
	e.erase(unique(e.begin(), e.end()), e.end());
}

void define_mapping(std::string& s, std::map<char, int>& m)
{
	for (size_t i = 0; i < s.size(); ++i)
		m.insert(std::pair<char, int>(s[i], i));
}

void define_mapping_fibre(std::vector<std::string>& s,
		std::map<std::string, int>& m)
{
	for (size_t i = 0; i < s.size(); ++i)
		m.insert(std::pair<std::string, int>(s[i], i));
}

void print_alignment(const Alignment<int>& alignment, const std::string& sigma1,
		const std::string& sigma2, const matching_schema<bool>& m,
		int& distance, bool is_identity)
{
	// print first string
	for (size_t i = 0; alignment.a[i] != _END_ALIGNMENT; i++)
		std::cout
				<< (alignment.a[i] != _GAP_FLAG ? sigma1[alignment.a[i]] : '-');
	std::cout << endl;

	// print second string
	for (size_t i = 0; alignment.b[i] != _END_ALIGNMENT; i++)
		std::cout
				<< (alignment.b[i] != _GAP_FLAG ? sigma2[alignment.b[i]] : '-');
	std::cout << endl;

	// print the alignment (a string containing * or space)
	size_t index = 0;
	size_t temp = 0;
	for (; alignment.a[index] != _END_ALIGNMENT; index++)
	{
		// if both chars are not gaps
		if (alignment.a[index] != _GAP_FLAG && alignment.b[index] != _GAP_FLAG)
		{
			if (!m.ms[alignment.a[index]][alignment.b[index]] || // if they are in match or
					(is_identity
							&& (sigma1[alignment.a[index]]
									== sigma2[alignment.b[index]])))
			{	// the self_identity is active and they are the same char
				std::cout << '*';
				temp++;
			}
			else
				std::cout << ' ';
		}
		else
			std::cout << ' ';
	}
	std::cout << endl;

	// here we return the new edit distance including the self_identity case
	distance = index - temp;
}
