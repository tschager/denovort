// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Thomas Tschager $
// $Authors: Yves Frank, Simon Rösch, Thomas Tschager, and Valentin Venzin $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_DENOVO_SYMDIFF_AMINOACIDS_H
#define OPENMS_ANALYSIS_DENOVO_SYMDIFF_AMINOACIDS_H

#include <string>
#include <map>
#include <cmath>
#include <set>
#include <fstream>
#include "data_types.h"
#include "config.h"

// Set containing pairs of single amino acids and all combinations of two amino acids and their masses
std::set<std::pair<double, std::string> > amino_masses;
// Same as @em amino_masses but 'I' and 'L' are considered separately
std::set<std::pair<double, std::string> > amino_masses_withI;
std::set<peak, comp_peak> spectrum;

// Array containing pairs single amino acids and their corresponding masses
std::pair<double, std::string> aminoacids[19] = 
{
  std::mp(57.02147,"G"), // Note that @em mp is a shortcut for @em make_pair
  std::mp(71.03712,"A"), 
  std::mp(87.03203,"S"), 
  std::mp(97.05277,"P"), 
  std::mp(99.06842,"V"), 
  std::mp(101.04768,"T"), 
  std::mp(103.00919,"C"), 
  std::mp(113.08407,"L"), 
  std::mp(114.04293,"N"), 
  std::mp(115.02695,"D"), 
  std::mp(128.05858,"Q"), 
  std::mp(128.09497,"K"), 
  std::mp(129.04260,"E"), 
  std::mp(131.04049,"M"), 
  std::mp(137.05891,"H"), 
  std::mp(147.06842,"F"), 
  std::mp(156.10112,"R"), 
  std::mp(163.06333,"Y"), 
  std::mp(186.07932,"W")
};

// Same as @em aminoacids but 'I' and 'L' are considered separately
std::pair<double, std::string> aminoacids_withI[20] = 
{
  std::mp(57.02147,"G"), 
  std::mp(71.03712,"A"), 
  std::mp(87.03203,"S"), 
  std::mp(97.05277,"P"), 
  std::mp(99.06842,"V"), 
  std::mp(101.04768,"T"), 
  std::mp(103.00919,"C"), 
  std::mp(113.08407,"L"), 
  std::mp(113.08407,"I"), 
  std::mp(114.04293,"N"), 
  std::mp(115.02695,"D"), 
  std::mp(128.05858,"Q"), 
  std::mp(128.09497,"K"), 
  std::mp(129.04260,"E"), 
  std::mp(131.04049,"M"), 
  std::mp(137.05891,"H"), 
  std::mp(147.06842,"F"), 
  std::mp(156.10112,"R"), 
  std::mp(163.06333,"Y"), 
  std::mp(186.07932,"W")
};

/**
  @return a lookup table @em AA_int2char (map) that maps a integers to a amino acids
*/
std::map<int, char> create_AA_int2char()
{
  std::map<int, char> AA_int2char;
  AA_int2char[ 0] = '-'; 
  AA_int2char[ 1] = 'G'; 
  AA_int2char[ 2] = 'A'; 
  AA_int2char[ 3] = 'S'; 
  AA_int2char[ 4] = 'P'; 
  AA_int2char[ 5] = 'V'; 
  AA_int2char[ 6] = 'T'; 
  AA_int2char[ 7] = 'C'; 
  AA_int2char[ 8] = 'L'; 
  AA_int2char[ 9] = 'I'; 
  AA_int2char[10] = 'N'; 
  AA_int2char[11] = 'D'; 
  AA_int2char[12] = 'Q'; 
  AA_int2char[13] = 'K'; 
  AA_int2char[14] = 'E'; 
  AA_int2char[15] = 'M'; 
  AA_int2char[16] = 'H'; 
  AA_int2char[17] = 'F'; 
  AA_int2char[18] = 'R'; 
  AA_int2char[19] = 'Y'; 
  AA_int2char[20] = 'W';
  return AA_int2char;
}
std::map<int, char> AA_int2char = create_AA_int2char();

/**
  @return a lookup table @em AA_char2int (map) that maps a amino acids to an integers
*/
std::map<char, int> create_AA_char2int()
{
  std::map<char, int> AA_char2int;
  AA_char2int['-'] =  0;
  AA_char2int['G'] =  1;
  AA_char2int['A'] =  2;
  AA_char2int['S'] =  3;
  AA_char2int['P'] =  4;
  AA_char2int['V'] =  5;
  AA_char2int['T'] =  6;
  AA_char2int['C'] =  7;
  AA_char2int['L'] =  8;
  AA_char2int['I'] =  9;
  AA_char2int['N'] = 10;
  AA_char2int['D'] = 11;
  AA_char2int['Q'] = 12;
  AA_char2int['K'] = 13;
  AA_char2int['E'] = 14;
  AA_char2int['M'] = 15;
  AA_char2int['H'] = 16;
  AA_char2int['F'] = 17;
  AA_char2int['R'] = 18;
  AA_char2int['Y'] = 19;
  AA_char2int['W'] = 20;
  return AA_char2int;
}
std::map<char, int> AA_char2int = create_AA_char2int();

/**
  @brief Initializes the tables @em amino_masses and @em amino_masses_withI. We insert masses of single amino
  acids and masses of pairs of amino acids to permit gaps in the list of prefixes.
  
  The tables @em amino_masses and @em amino_masses_withI can be used to find the labels for the edges of the 
  spectrum graph.
*/
void init_mass_table() 
{
  if (amino_masses.size() == 0)
  {
    for (int i = 0; i < 19; ++i)
    {
      amino_masses.insert(mp(aminoacids[i].first, aminoacids[i].second));
    }
    for (int i = 0; i < 19; ++i) 
    {
      for (int j = 0; j < 19; ++j) 
      {
        std::pair<double,std::string> p = std::mp(aminoacids[i].first + aminoacids[j].first,  "("+aminoacids[i].second + aminoacids[j].second+")");
        amino_masses.insert(p);
      }
    }
  }

  if (amino_masses_withI.size() == 0)
  {
    for (int i = 0; i < 20; ++i)
    {
      amino_masses_withI.insert(mp(aminoacids_withI[i].first, aminoacids_withI[i].second));
    }
    for (int i = 0; i < 20; ++i) 
    {
      for (int j = 0; j < 20; ++j) 
      {
        std::pair<double,std::string> p = std::mp(aminoacids_withI[i].first + aminoacids_withI[j].first, "(" + aminoacids_withI[i].second + aminoacids_withI[j].second + ")");
        amino_masses_withI.insert(p);
      }
    }
  }
}

/**
  @brief Computes list of strings (of amino acids) that match a given mass. This corresponds to finding
  the edge labels of the spectrum graph.

  @param mass Mass difference between two prefix masses - two peaks in the spectrum
  @param conf Parameters for algorithm
  @return list @em res of strings matching the mass
*/
std::vector<std::string> get_label(double mass, config conf) 
{
  std::vector<std::string> res;
  std::set<std::pair<double, std::string> >& masses = (conf.ENABLE_I) ? amino_masses_withI : amino_masses;

  std::set<std::pair<double, std::string> >::iterator it = masses.lower_bound(std::mp(mass - conf.EPS, ""));
  for (; it != masses.end() && (*it).first - mass <= conf.EPS; ++it) 
  {
    if (fabs((*it).first - mass) <= conf.EPS) 
      res.push_back((*it).second);
  }
  return res;
}

/**
  @return the mass of the amino acid @em c
*/
double get_aminoacid_mass(char c)
{
  for (int i = 0; i < sizeof(aminoacids_withI) / sizeof(aminoacids_withI[0]); i++)
  {
    std::pair<double, std::string> p = aminoacids_withI[i];
    if (p.second[0] == c)
      return p.first;
  }
  return 0;
}

/**
  @brief Get all prefix masses for a sequence of amino acids.

  @param seq Sequence of amino acids
  @return list of prefix masses
*/
std::vector<double> get_prefix_masses(std::string seq) 
{
  std::vector<double> prefix_masses;
  for (int i = 0; i < seq.size(); i++)
  {
    std::string aa = seq.substr(i, 1);
    // substitute all L's by I's
    if (aa == "L")
      aa == "I";
    if (prefix_masses.size() == 0)
      prefix_masses.push_back(get_aminoacid_mass(aa.at(0)));
    else
      prefix_masses.push_back(prefix_masses.back() + get_aminoacid_mass(aa.at(0)));
  }
  return prefix_masses;
}

/**
  @brief Computes the recall of an amino acid sequence @em seq2 with respect to the true amino acid sequence 
  @em seq1.

  @param seq1 True amino acid sequence
  @param seq2 Some other amino acid sequence
  @param eps Tolerance parameter (>0)
  @return recall = number of common prefix masses of @em seq1 and @em seq2 divided by number of prefix masses 
  of @emseq1
*/  
double get_recall(std::string seq1, std::string seq2, const double eps)
{
  std::vector<double> prefix_masses_seq1 = get_prefix_masses(seq1);
  std::vector<double> prefix_masses_seq2 = get_prefix_masses(seq2);
  std::set<double> intersection(prefix_masses_seq1.begin(), prefix_masses_seq1.end());

  bool found = true;
  double low,high;
  while (found)
  {
    found = false;
    std::set<double>::iterator iter;
    for (iter = intersection.begin(); iter != intersection.end(); ++iter) 
    {
      double m = *iter;
      std::vector<double>::iterator up = upper_bound(prefix_masses_seq2.begin(), prefix_masses_seq2.end(), m);
      if (up != prefix_masses_seq2.begin()) 
        low = *(up - 1);
      else
        low = prefix_masses_seq2.front();
      high = *up;
      if (fabs(high - m) > eps && fabs(low - m) > eps)
      {
        intersection.erase(m);
        found = true;
      }
    }
  }

  double recall = (double)intersection.size() / prefix_masses_seq1.size();
  return recall;
}

#endif