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

#ifndef OPENMS_ANALYSIS_DENOVO_SYMDIFF_SOLVEWRAPPER_H
#define OPENMS_ANALYSIS_DENOVO_SYMDIFF_SOLVEWRAPPER_H

#include <string>
#include <iostream>
#include <map>
#include <algorithm>
#include "config.h"

#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace OpenMS;

#include "symdiffscore_multyions.h"
#include "amino_acids.h"
#include "peptide_solver_symmetric_difference_general.h"

#include "retention_time.h"
#include "PeptideHitRT.h"
const double H_MASS = 1.00794; // The mass of H molecule in Daltons
const double H20_MASS = 18.01048; // The mass of H20 molecule in Daltons

// Used to sort the peptides by their score
bool peptideOrder(peptide a, peptide b) 
{ 
  return a.first > b.first || (a.first == b.first && a.second > b.second); 
}

/**
  @brief Struct that represents a peak pair and is needed in the merge process
  @Note We need @em low_mass_h and @em high_mass_h to reconstruct spectrum in full range (0 to W) given a merged
  spectrum only in range 0 to W/2.
*/
typedef struct{
  double time;
  double mz;
  double h;
  double low_mass_h;
  double high_mass_h;
  bool deleted;
} mergepeak;


/**
  @brief Merge similar (close) peaks in the mass spectrum. After preprocessing, search routine is called.

  @param spectrum List of peaks in the MS
  @param mass of the peptide in question
  @param time Measured retention time of peptide in question
  @param offset Mass offset of the last character
  @param result A list of peptides to be populated by the search routine
  @param last_letter (Given) last char of peptide in question
  @param rho Relative score a solution must have to be listed in the results
  @param conf Parameters for algorithm
  @param csv_path Path to a csv file containing the RT coefficients that are to be used
*/
void solve_prec_mass_with_offset(std::vector<peak>& spectrum, double mass, double time, double offset, std::vector<peptide>& result, char last_letter, double rho, config conf, std::string csv_path)
{

  if (spectrum.size() == 0) 
    return;

  // Mass of the amino acid sequence without the C-terminus and without the tag of the last amino acid
  double W = mass - offset - H20_MASS - H_MASS;
  RetentionTime rt = RetentionTime(csv_path, W, time, conf);

  // Throw away masses that are too 'heavy'
  std::vector<peak> v_spectrum;
  for (int i = 0; i < spectrum.size(); i++)
  {
    if (spectrum[i].mz < mass)
      v_spectrum.push_back(spectrum[i]);
  }

  /**
    Generate all potential b-ion-masses (uncharged) in range 0 - (mass - H_MASS)/2. 
  */

  std::vector<mergepeak> potential_b_peaks;
  // Contains pairs (intensity, index) s.t. the peaks are sorted wrt. intensity
  std::set<std::pair<double, int> > peakIntensityOrdering; 
  int front = 0;
  int back = v_spectrum.size()-1;

  while (front <= back) 
  {
    if (mass - v_spectrum[back].mz > v_spectrum[front].mz - H_MASS)
    {
      peak p = v_spectrum[front];
      mergepeak mpeak;
      mpeak.time = p.time;
      mpeak.mz = p.mz - H_MASS;
      mpeak.h = p.h;
      mpeak.low_mass_h = p.h;
      mpeak.high_mass_h = 0;
      mpeak.deleted = false;
      peakIntensityOrdering.insert(std::make_pair(-p.h, potential_b_peaks.size()));
      potential_b_peaks.push_back(mpeak);
      front++;
    }
    else
    {
      peak p = v_spectrum[back];
      mergepeak mpeak;
      mpeak.time = p.time;
      mpeak.mz = mass - p.mz;
      mpeak.h = p.h;
      mpeak.low_mass_h = 0;
      mpeak.high_mass_h = p.h;
      mpeak.deleted = false;
      peakIntensityOrdering.insert(std::make_pair(-p.h, potential_b_peaks.size()));
      potential_b_peaks.push_back(mpeak);
      back--;
    }
  }

  /** 
    Greedy algorithm for peak merging: Go through peaks in decreasing intensity order. For every peak merge
    with closest neighbor and update position as weighted average, repeat this until no two peaks are closer
    than the merging accuracy.
  */
  for (std::set<std::pair<double, int> >::iterator iter = peakIntensityOrdering.begin(); iter != peakIntensityOrdering.end(); iter++)
  {
    std::pair<double, int> peakPos = *iter;
    int pos = peakPos.second; // Index of current peak
    if (potential_b_peaks[pos].deleted)
      continue;

    // Indices of neighboring peaks of current peak
    int left = pos - 1;
    int right = pos + 1;
    
    // Merge current peak with undeleted neighbor, if close enough and not deleted
    while ((left >= 0 && fabs(potential_b_peaks[left].mz - potential_b_peaks[pos].mz) < conf.MERGING_EPS) || (right < potential_b_peaks.size() && fabs(potential_b_peaks[right].mz - potential_b_peaks[pos].mz) < conf.MERGING_EPS))
    {
      // Compute mass differences to left and right neighbor
      double leftDiff = conf.MERGING_EPS;
      if (left >= 0)
        leftDiff = fabs(potential_b_peaks[left].mz - potential_b_peaks[pos].mz);
      double rightDiff = conf.MERGING_EPS;
      if (right < potential_b_peaks.size())
        rightDiff = fabs(potential_b_peaks[right].mz - potential_b_peaks[pos].mz);

      if (leftDiff < rightDiff)
      {
        // Merge current with left neighbor if not deleted
        if (!potential_b_peaks[left].deleted)
        {
          potential_b_peaks[pos].mz = (potential_b_peaks[pos].mz * potential_b_peaks[pos].h + potential_b_peaks[left].mz * potential_b_peaks[left].h) / (potential_b_peaks[pos].h + potential_b_peaks[left].h);
          potential_b_peaks[pos].h += potential_b_peaks[left].h;
          potential_b_peaks[pos].low_mass_h += potential_b_peaks[left].low_mass_h;
          potential_b_peaks[pos].high_mass_h += potential_b_peaks[left].high_mass_h;
          potential_b_peaks[left].deleted = true;
        }
        left--;
      }
      else
      {
        // Merge current with right neighbor if not deleted
        if (!potential_b_peaks[right].deleted)
        {
          potential_b_peaks[pos].mz = (potential_b_peaks[pos].mz * potential_b_peaks[pos].h + potential_b_peaks[right].mz * potential_b_peaks[right].h) / (potential_b_peaks[pos].h + potential_b_peaks[right].h);
          potential_b_peaks[pos].h += potential_b_peaks[right].h;
          potential_b_peaks[pos].low_mass_h += potential_b_peaks[right].low_mass_h;
          potential_b_peaks[pos].high_mass_h += potential_b_peaks[right].high_mass_h;
          potential_b_peaks[right].deleted = true;
        }
        right++;
      }
    }
  }

  /**
    Create list of peaks in whole range and only keep peaks with h > NOISE_H_CUTOFF. This means that each
    merged peak in potential_b_peaks gets added twice to @em peaks; once on the 'left side' of the spectrum
    and once on the 'right side'. Restart this procedure and double the NOISE_H_CUTOFF until we get a spectrum
    with less than 500 peaks.
  */

  std::vector<peak> peaks;
  double cutoff = 0;

  do 
  {
    peaks.clear();
    peak startpeak = {0, 0, 0};
    peaks.push_back(startpeak);

    for (int i = 0; i < potential_b_peaks.size(); i++)
    {
      mergepeak p = potential_b_peaks[i];
      if (!p.deleted && p.h > cutoff)
      {
        peak pp;
        pp.time = p.time;
        pp.mz = p.mz;
        pp.h = p.low_mass_h;
        peaks.push_back(pp);
      }
    }

    for (int i = potential_b_peaks.size() - 1; i >= 0; i--)
    {
      mergepeak p = potential_b_peaks[i];
      if (!p.deleted && p.h > cutoff)
      {
        peak pp;
        pp.time = p.time;
        pp.mz = (mass - p.mz) - 1 * H_MASS;
        pp.h = p.high_mass_h;
        peaks.push_back(pp);
      }
    }

    peak endpeak = {0, W, 0};
    peaks.push_back(endpeak);

    cutoff += conf.NOISE_H_CUTOFF;

  } while (peaks.size() > 500);

  // Call the search routine that finds amino acid sequences with preprocessed spectrum
  find_all_peptides(peaks, last_letter, rt, rho, result, conf);
}

/**
  @brief Given the spectrometry data, this method invokes the search procedure for every possible last 
  character. Afterwards all solutions are added to a PeptideIdentification (list of possible peptide hits for
  this spectrum) which is then returned.

  @param spectrum List of peaks in the MS
  @param mass of the peptide in question
  @param rho Relative score a solution must have to be listed in the results
  @param time Measured retention time of the peptide in question
  @param conf Parameters for algorithm
  @param lastChars The last character of a peptide can only be one of these characters
  @param lastCharTags Mass of amino acids in @em lastChars
  @param csv_path Path to a csv file containing the RT coefficients that are to be used
  @return All peptide hits that are at least @em rho times as good as the highest scoring hit
*/
PeptideIdentification solve_one_prec_mass(std::vector<peak>& spectrum, double mass, double rho, double time, config conf, StringList lastChars, DoubleList lastCharTags, std::string csv_path) 
{

  init_mass_table();
  std::vector<peptide> peptides; // List of results; needs to be post-processed
  
  // Call search procedure for for each possible last character
  for (int i = 0; i < lastChars.size(); i++)
  {
    solve_prec_mass_with_offset(spectrum, mass, time, lastCharTags[i], peptides, lastChars[i][0], rho, conf, csv_path);
  }

  /**
    Sort, prune and bring output, i.e. @em peptides, into framework compatible form
  */

  PeptideIdentification identification;
  identification.setRT(time);
  identification.setHigherScoreBetter(true);
  identification.setSignificanceThreshold(rho);

  if (peptides.size() > 0)
  {
    int max_score = 0;
    std::sort(peptides.begin(), peptides.end(), peptideOrder); // Sort descending
    max_score = std::max(max_score, peptides[0].first.first);

    for (int j = 0; j < peptides.size(); j++)
    {
      peptide p = peptides[j];
      // Only keep peptides that have score at least @em rho * max_score - Note that this step is necessary
      // because we have results from searches for multiple last characters
      if (p.first.first >= rho * max_score)
      {
        // Generate peptide hit (remove brackets of multi amino acid edges)
        std::stringstream ss;
        for (int si = 0; si < p.second.size(); si++)
        {
          for (int ci = 0; ci < p.second[si].size(); ci++)
          {
            char c = p.second[si][ci];
            if (c != '(' && c != ')')
            {
              ss << c;
            }
          }
        }
        PeptideHitRT hit = PeptideHitRT(p.first.first, 0, 1, AASequence::fromString(ss.str()), p.first.second);
        identification.insertHit(hit);
      }
    }
    identification.assignRanks();
  }
  return identification;
}

#endif
