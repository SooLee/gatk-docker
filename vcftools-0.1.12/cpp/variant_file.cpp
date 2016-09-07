/*
 * variant_file.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: amarcketta
 */

#include "variant_file.h"

variant_file::~variant_file() {}

// Return the number of individuals that have not been filtered out
int variant_file::N_kept_individuals() const
{
	int N_kept = 0;
	for (unsigned int ui=0; ui<include_indv.size(); ui++)
		if (include_indv[ui] == true)
			N_kept++;
	return N_kept;
}

// Return the number of sites that have not been filtered out
int variant_file::N_kept_sites() const
{
	return N_kept_entries;
}

// Return the total number of sites in the file
int variant_file::N_total_sites() const
{
	return N_entries;
}

void variant_file::ByteSwap(unsigned char *b, int n) const
{
   register int i = 0;
   register int j = n-1;
   while (i<j)
   {
      std::swap(b[i], b[j]);
      i++, j--;
   }
}

void variant_file::get_default_contigs(vector<string> &contig_vector)
{
	contig_vector.resize(0);
	contig_vector.push_back("##contig=<ID=1,length=249250621,assembly=b37>");
	contig_vector.push_back("##contig=<ID=2,length=243199373,assembly=b37>");
	contig_vector.push_back("##contig=<ID=3,length=198022430,assembly=b37>");
	contig_vector.push_back("##contig=<ID=4,length=191154276,assembly=b37>");
	contig_vector.push_back("##contig=<ID=5,length=180915260,assembly=b37>");
	contig_vector.push_back("##contig=<ID=6,length=171115067,assembly=b37>");
	contig_vector.push_back("##contig=<ID=7,length=159138663,assembly=b37>");
	contig_vector.push_back("##contig=<ID=8,length=146364022,assembly=b37>");
	contig_vector.push_back("##contig=<ID=9,length=141213431,assembly=b37>");
	contig_vector.push_back("##contig=<ID=10,length=135534747,assembly=b37>");
	contig_vector.push_back("##contig=<ID=11,length=135006516,assembly=b37>");
	contig_vector.push_back("##contig=<ID=12,length=133851895,assembly=b37>");
	contig_vector.push_back("##contig=<ID=13,length=115169878,assembly=b37>");
	contig_vector.push_back("##contig=<ID=14,length=107349540,assembly=b37>");
	contig_vector.push_back("##contig=<ID=15,length=102531392,assembly=b37>");
	contig_vector.push_back("##contig=<ID=16,length=90354753,assembly=b37>");
	contig_vector.push_back("##contig=<ID=17,length=81195210,assembly=b37>");
	contig_vector.push_back("##contig=<ID=18,length=78077248,assembly=b37>");
	contig_vector.push_back("##contig=<ID=19,length=59128983,assembly=b37>");
	contig_vector.push_back("##contig=<ID=20,length=63025520,assembly=b37>");
	contig_vector.push_back("##contig=<ID=21,length=48129895,assembly=b37>");
	contig_vector.push_back("##contig=<ID=22,length=51304566,assembly=b37>");
	contig_vector.push_back("##contig=<ID=X,length=155270560,assembly=b37>");
	contig_vector.push_back("##contig=<ID=Y,length=59373566,assembly=b37>");
	contig_vector.push_back("##contig=<ID=MT,length=16569,assembly=b37>");
}
