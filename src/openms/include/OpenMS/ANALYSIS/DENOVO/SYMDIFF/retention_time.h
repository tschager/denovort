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

#ifndef OPENMS_ANALYSIS_DENOVO_SYMDIFF_RETENTIONTIME_H
#define OPENMS_ANALYSIS_DENOVO_SYMDIFF_RETENTIONTIME_H

#include <string>
#include <map>
#include <cmath>
#include <fstream>
#include "data_types.h"
#include "config.h"
#include "amino_acids.h"

#include <limits>
#include <cfloat>

int gamma_l; // Number of different position dependent RT coefficients for the prefix
int gamma_r; // Number of different position dependent RT coefficients for the suffix
bool model_lin;
bool model_nei;
bool model_pos;
// Map that stores all coefficients (depending on the model) of each amino acid
std::map<char, std::vector<int> > coefs;

/**
  @brief Removes brackets "(" and ")" at beginning and end of @em seq.

  @param seq Sequence to be stripped of brackets
  @return sequence without brackets
*/
std::string removeBrackets(std::string seq)
{
  if (seq[0] == '(' && seq[seq.size() - 1] == ')') 
    seq = seq.substr(1, seq.size() - 2);
  return seq;
}

/**
  @brief Computes retention time of an amino acid sequence @em seq with given positions of the leftmost and
  the rightmost characters.

  @Note Must be called after initialization of coefs

  @param seq Amino acid sequence
  @param pos_l Position of the leftmost character of @em seq in the peptide. If position is unknown, set to
    some value greater than @em gamma_l
  @param pos_r Position of the rightmost character of @em seq in the peptide. If position is unknown, set to
    some value greater than @em gamma_r
  @return predicted retention time of @em seq
*/
int get_seq_rt(std::string seq, int pos_l, int pos_r)
{
  // Remove brackets of multi-char labels
  if (seq[0] == '(' && seq[seq.size() - 1] == ')') 
    seq = seq.substr(1, seq.size() - 2);
  
  int res = 0;
  int n = seq.length();
  int c_l = std::min(pos_l, gamma_l);
  int c_r = std::min(pos_r, gamma_r);
  int l = 0; int r = n - 1;
  
  // Left position dependent coefficient
  while (c_l < gamma_l && l <= r) 
  {
    int plus = coefs[seq[l]][1+c_l];
    res += plus;
    l++; c_l++;	
  }
  // Right position dependent coefficient
  while (c_r < gamma_r && l <= r) 
  {
    int plus = coefs[seq[r]][gamma_l + gamma_r - c_r];
    res += plus;
    r--; c_r++;	
  }
  // Default coefficients
  for (int i = l; i <= r; i++) 
  {
    res += coefs[seq[l++]][0];
  }
  return res;
}

/**
  Computes retention time of seq in the neighborhood rt prediction model for a sequence seq.

  Must be called after initialization of coefs
*/
int get_seq_rt_nei(std::string seq)
{
  // Remove brackets of multi-char labels
  if (seq[0] == '(' && seq[seq.size() - 1] == ')') 
    seq = seq.substr(1,seq.size() - 2);

  int res = 0;
  int n = seq.length();
  for (int i = 0; i < n - 1; i++) 
  {
    res += coefs[seq[i]].at(AA_char2int.at(seq[i + 1]));
  }
  return res;
}

/**
	@brief Initialize the retention time coefficients. This fills map @em coefs. If it is already initialized
  nothing happens.
	
  @param csv_path Path to csv file containing RT coefficients in the correct format (minimal error handling only)

  @exception Exception::InvalidParameter thrown if the path to the csv file is invalid, the first line is not
  @em "min" or @em "sec", the second line is not @em "lin", @em "pos", or @em "nei"
*/
void init_coefs(std::string csv_path)
{
  if (coefs.empty()) 
  {
    std::ifstream coef_file;
    coef_file.open(csv_path.c_str());
    if (!coef_file.good()) 
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Invalid path to .csv file.");
    
    std::string line;
    getline(coef_file, line, ',');
    if (line != "sec" && line != "min") 
    {
      std::string err = "csv file has incorrect form; expected 'sec' or 'min', got " + line;
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, err);
    }
    bool minutes = false;
    if (line == "min")
      minutes = true;
    
    model_lin = false;
    model_pos = false;
    model_nei = false;
    
    getline(coef_file, line, ',');
    if (line != "lin" && line != "pos" && line != "nei")
    {
      std::string err = "wrong csv file format; model declaration missing; expected 'lin', 'pos' or 'nei', got " + line;
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, err);
    }
    
    // Read position dependent coefficients
    if (line == "pos")
    {
      LOG("Position-dependent RT prediction model \n");
      model_pos = true;
      getline(coef_file, line, ',');
      getline(coef_file, line, ','); // gamma_l value
      std::istringstream(line) >> gamma_l;	
      getline(coef_file, line, ',');
      getline(coef_file, line, ','); // gamma_r value
      std::istringstream(line) >> gamma_r;
      
      // Skip header
      for (int i = 0; i < 2 + gamma_l + gamma_r; i++) 
      {
        getline(coef_file, line, ',');
      }

      while (coef_file.good()) 
      {
        // Contains @em gamma_l + @em gamma_r position dependent RT coefficients and one default coefficient
        std::vector<int> aa_coefs(1 + gamma_l + gamma_r, 0.0);
        std::string aa; // 
        getline(coef_file, aa, ',');
        if (aa[1] == '\0')
          break;
        getline(coef_file, line, ','); // Default
        std::istringstream(line) >> aa_coefs[0];
        for (int i = 0; i < gamma_l; i++) // @em gamma_l prefix coefficients
        {
          getline(coef_file, line, ',');
          std::istringstream(line) >> aa_coefs[1 + i];
        }
        for (int i = 0; i < gamma_r; i++) // @em gamma_r suffix coefficients
        {
          getline(coef_file, line, ',');	
          std::istringstream(line) >> aa_coefs[1 + gamma_l + i];		
        }

        // Convert Coefficients from minutes to seconds if necessary
        if (minutes)
        {
          for (int i = 0; i < aa_coefs.size(); i++) 
          {
            aa_coefs[i] *= 60;
          }
        }
        coefs[aa[1]] = aa_coefs;
      }
    }
    // Read neighboring coefficients
    else if (line == "nei")
    {
      LOG("Neighborhood-dependent RT prediction model \n");
      model_nei = true;
      gamma_l = 0;
      gamma_r = 0;
      
      // Skip header
      for (int i = 0; i < 3; i++) 
      {
        getline(coef_file, line, ',');
      }

      while (coef_file.good()) 
      {
        std::string aa;
        std::string aa_prev;
        getline(coef_file, aa, ',');
        if (aa[1] == '\0')
          break;
        getline(coef_file, aa_prev, ',');
        
        if (AA_char2int.find(aa_prev[0]) == AA_char2int.end()) 
        {
          std::string err = "wrong csv file format; previous amino acid character unknown, got " + std::string(1, aa_prev[1]);
          throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, err);
        }
        else 
        {
          getline(coef_file, line, ',');
          
          if (coefs.find(aa[1]) == coefs.end())
          {
            // First coefficient of amino acid aa[1]
            std::vector<int> aa_coefs(21, 0.0);
            std::istringstream(line) >> aa_coefs[AA_char2int.at(aa_prev[0])];
            coefs[aa[1]] = aa_coefs;
          }
          else
            std::istringstream(line) >> coefs[aa[1]][AA_char2int.at(aa_prev[0])];
        }

        // Convert Coefficients from minutes to seconds if necessary
        if (minutes)
        {
          for (int i = 0; i < coefs[aa[1]].size(); i++) 
          {
            coefs[aa[1]][i] *= 60;
          }
        }
      }
    }
    // Read linear coefficients
    else
    {
      LOG("Linear RT prediction model \n");
      model_lin = true;
      gamma_l = 0;
      gamma_r = 0;
      
      // Skip header
      for (int i = 0; i < 2; i++) 
      {
        getline(coef_file, line, ',');
      }

      while (coef_file.good()) 
      {
        std::vector<int> aa_coefs(1, 0.0);
        std::string aa;
        std::string key;
        getline(coef_file, aa, ',');
        if (aa[2] == '\0')
          break;
        getline(coef_file, line, ','); // Default
        std::istringstream(line) >> aa_coefs[0];
        
        // Convert Coefficients from minutes to seconds
        if (minutes)
        {
          for (int i = 0; i < aa_coefs.size(); i++) 
          {
            aa_coefs[i] *= 60;
          }
        }
        coefs[aa[2]] = aa_coefs;
      }
    }
    coef_file.close();
  }
}

/**
	@brief Class that models the retention time of a peptide. The exact retention time is rounded to the nearest
  integer. Furthermore the smallest and the largest possible retention time for a peptide with a given mass
  is calculated.

  Use @em get_seq_rt() and @em get_seq_rt_nei() to compute the retention time of a specific sequence rather
  than this class.
*/	
class RetentionTime {
  config cnf;
  double exact_time; 
  int time, t_min, t_max;	
  double mass; // Mass of sequence without the C-terminus and without the tag of the last amino acid
  public:
    double get_exact_rt() { return exact_time; }
    int get_rt() { return time; }
    int get_t_min() { return t_min; } // Get smallest possible retention time for peptide with mass @em mass
    int get_t_max() { return t_max; } // Get largest possible retention time for peptide with mass @em mass
    std::pair<int,int> get_gamma() { return std::mp(gamma_l,gamma_r); }
    bool isLinModel() { return model_lin; }
    bool isPosModel() { return model_pos; }
    bool isNeiModel() { return model_nei; }
    RetentionTime(std::string csv_path, double W, double T, config conf);
    RetentionTime();
};

/**
  @brief Constructor.

  @param csv_path Path to file containing retention time coefficients
  @param W Measured time of the peptide
  @oaram config Parameters for the algorithm
*/
RetentionTime::RetentionTime(std::string csv_path, double W, double T, config conf)
{
  cnf = conf;
  exact_time = T;
  time = round(exact_time);
  mass = W;
  
  // Read RT coefficients from csv file
  LOG("Read RT coefficient file " + csv_path + "\n");
  init_coefs(csv_path);

  /**
    Compute smallest and largest possible retention time for peptide with mass @em mass. 
    A peptide would have retention time @em t_min if it would only consist of the amino acid with the smallest
    retention time to mass ratio. Analogous for @t_max.
  */
  double min_frac = DBL_MAX;
  double max_frac = DBL_MIN;
  LOG("Begin\n");
  for (auto a : AA_char2int) 
  {
    if (coefs.find(a.first) == coefs.end())
    {
      LOG(" *** WARNING no coefs for " + std::string(1, a.first) + "\n");
    }
    else
    {
      if (coefs.at(a.first).size() == 0 || a.first == '-') //we ignore the character '-' (begin and end of sequence)
        continue;
      double max = *max_element(std::begin(coefs.at(a.first)), std::end(coefs.at(a.first)));
      double min = *min_element(std::begin(coefs.at(a.first)), std::end(coefs.at(a.first)));
      min_frac = std::min(min_frac, min / get_aminoacid_mass(a.first));
      max_frac = std::max(max_frac, max / get_aminoacid_mass(a.first));
    }
  }
  t_min = floor(min_frac * mass);
  t_max = ceil(max_frac * mass);
}

RetentionTime::RetentionTime() {}

/**
	@brief This class represents all retention time values that may be added to a given retention time value
	such that their sum is a feasible retention time. I.e. given t, all t' with |t+t'-T|<RT_EPS are feasible.

  This is useful if we want to extend a prefix of a peptide. For the extension of the prefix, we only need
  to consider retention time values between @em get_min_feasible() and @em get_max_feasible().
*/
class FeasibleTimes {
  int t_target, ft_min, ft_max;
  public:
    FeasibleTimes(int t, RetentionTime rt, config conf) 
    {
      int meas_time = rt.get_rt(); // Rounded retention time of peptide
      t_target = meas_time - t;
      int tmp = t_target - conf.RT_EPS / 2;
      ft_min = std::max(tmp, rt.get_t_min()); // Round and check that does not exceed t_min
      if (ft_min < tmp)
        ft_min++;	
      tmp = t_target + conf.RT_EPS / 2;
      ft_max = std::min(tmp, rt.get_t_max());
      if (ft_max > tmp)
        ft_max--;
      
      LOG("Target retention time is " + std::to_string(rt.get_rt()) + " (t_target: " + std::to_string(t_target) + ")\n");
      LOG("RT_EPS is " + std::to_string(conf.RT_EPS / 2) + " and feasible retention times between " + std::to_string(ft_min) + " and " + std::to_string(ft_max) +"\n");
    }
    int get_min_feasible() { return ft_min; }
    int get_max_feasible() { return ft_max; }
    int get_nr_feasible_rts() {return ft_max - ft_min; }
    bool is_feasible(int t) { return (ft_min <= t && ft_max >= t); } // Check if @em t is in feasible RT interval
};

#endif
