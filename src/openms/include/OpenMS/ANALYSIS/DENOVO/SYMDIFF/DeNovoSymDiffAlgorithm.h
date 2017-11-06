// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#ifndef OPENMS_ANALYSIS_DENOVO_SYMDIFF_SYMDIFFDENOVOALGORITHM_H
#define OPENMS_ANALYSIS_DENOVO_SYMDIFF_SYMDIFFDENOVOALGORITHM_H

#define DEBUG
#ifdef DEBUG
#	define 	LOG(x) std::clog << x;
#else
#	define LOG(x) do {} while (0)
#endif

// OpenMS includes
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>

#include "solve_wrapper.h"
#include "config.h"
#include "peptide_solver_symmetric_difference_general.h"
#include "symdiffscore_multyions.h"
#include "retention_time.h"

namespace OpenMS
{
  /**
	 @brief A class handling the parameters to run the Retention Time Symmetric Difference DeNovo sequencing
    algorithm.

	 This class defines all parameters for the algorithm and it's default values. The parameters are filled into
    a config struct and passed to the solver. The only parameter check that is conducted, is that the number of 
    @em last_characters matches the number of @em last_character_tags. A basic parameter check is already done 
    by the framework. All other parameters are assumed to be correct.
  */

  class OPENMS_DLLAPI DeNovoSymDiffAlgorithm :
    public DefaultParamHandler
  {

public:

    /**
      @brief Default constructor that fills the default parameters
    */
    DeNovoSymDiffAlgorithm() :
      DefaultParamHandler("DeNovoSymDiffAlgorithm")
    {
      defaults_.setValue("rho", 0.9, "Relative score a solution must have to be listed in the results");

      StringList lastChars;
      lastChars.push_back("K");
      lastChars.push_back("R");
      defaults_.setValue("last_characters", lastChars, "The last character of a peptide can only be one of these characters");
      DoubleList lastCharTags;
      lastCharTags.push_back(8.0142);
      lastCharTags.push_back(10.008);
      defaults_.setValue("last_character_tag", lastCharTags, "For every last character a mass offset must be defined, use 0 if no mass offset is necessary");

      defaults_.setValue("accuracy", 0.02, "Accuracy used for equality of peaks (in Daltons)");
      defaults_.setValue("merge_accuracy", 0.04, "Accuracy used for merging peaks (in Daltons)");
      defaults_.setValue("noise_cutoff", 100.0, "Peaks with lower intensity are ignored");

      // Symdiff parameters
      defaults_.setValue("count_peaks", "false", "Count number of peaks instead of summing up their intensities (if true, missing_peak_punishment should be -2 and half_missing_peak_punishment should be -1 )");
      defaults_.setValidStrings("count_peaks", ListUtils::create<String>("true,false"));
      defaults_.setValue("half_missing_peak_punishment", -2500.0, "Punishment if only one peak of a b-/y-ion peak pair is missing");
      defaults_.setValue("missing_peak_punishment", -5000.0, "Negative score for a missing b-/y-ion peak pair");
      defaults_.setValue("missing_peak_threshold", 1.0, "Maximum intensity a peak can have to be counted as missing");
      defaults_.setValue("missing_ion_peak_punishment", -1000.0, "Punishment if a ion other than b-/y-ion is missing");

      // RT parameters
      defaults_.setValue("rt_accuracy", 1, "Discretization step in seconds taken for the retention time");
      defaults_.setValue("rt_window", 2500, "Total size (in seconds) of the time window around the measured retention time where the predicted retention time should be in order to qualify its peptide for a solution");
      defaults_.setValue("coef_path", "", "Directory where the .csv file with the retention time coefficients is located");
      defaults_.setValue("coef_file", "pos_dependent_rt.csv", "Name of the .csv file with the retention time coefficients");
      defaults_.setValue("enable_i", "false", "Boolean that determines whether L should be distinguished from I");

      // Multiple ion type parameters
      DoubleList b_ion_offsets;
      b_ion_offsets.push_back(-27.9949); // a-ion
      b_ion_offsets.push_back(-18.0106); // b-ion - H20
      defaults_.setValue("b_ion_offsets", b_ion_offsets, "Offsets for additional ion types to the b-ion");
      DoubleList y_ion_offsets;
      y_ion_offsets.push_back(-18.0106); // y-ion - H20
      y_ion_offsets.push_back(1); // y-ion with 13C isotope
      defaults_.setValue("y_ion_offsets", y_ion_offsets, "Offsets for additional ion types to the y-ion");

      defaults_.setValue("annotation_sequence", "", "A amino acid sequence from annotation data. This can be used for testing the algorithm with a test set where the correct sequence is already known. The tool then outputs performance values as well");
      
      // Write defaults into Param object param_
      defaultsToParam_();
    }

    /**
      @brief Execute a search for amino acid sequences given the mass spectrometer data. Output results.
      
      @param input Spectrum on which the search as to be performed
      @param mass Mass of peptide in question
      @Note Remaining parameters are passed by the framework over the @em param_ field.
      @result All peptide hits, provided their score is high enough, for supplied spectrum
      
      @exception Exception::InvalidParameter is thrown if the number of elements in the parameter 
      last_characters does not match the number of elements in last_character_tag.
    */
    PeptideIdentification searchSequences(const MSSpectrum<Peak1D>& input, double mass)
    {
      std::vector<peak> spectrum;

      for (Size i = 0; i < input.size(); i++)
      {
        peak p = {input.getRT(), input[i].getMZ(), input[i].getIntensity()};
        spectrum.push_back(p);
      }

      // Setup conf struct based on param_
      config conf;
      conf.EPS = param_.getValue("accuracy");
      conf.MERGING_EPS = param_.getValue("merge_accuracy");
      conf.NOISE_H_CUTOFF = param_.getValue("noise_cutoff");
      conf.COUNT_PEAKS = param_.getValue("count_peaks").toBool();
      conf.MISSING_PEAK_PUNISHMENT = param_.getValue("missing_peak_punishment");
      conf.HALF_MISSING_PEAK_PUNISHMENT = param_.getValue("half_missing_peak_punishment");
      conf.MISSING_PEAK_THRESHOLD_SCORE = param_.getValue("missing_peak_threshold");
      conf.MISSING_ION_PEAK_PUNISHMENT = param_.getValue("missing_ion_peak_punishment");
      conf.START_SCORE = -10 * conf.MISSING_PEAK_PUNISHMENT;
      conf.b_ion_offsets = param_.getValue("b_ion_offsets");
      conf.y_ion_offsets = param_.getValue("y_ion_offsets");

      StringList lastChars = param_.getValue("last_characters");
      DoubleList lastCharTags = param_.getValue("last_character_tag");

      if (lastChars.size()!= lastCharTags.size())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "For every last character a corresponding tag is required! Use 0 if you don't need an offset.");
      }

      conf.RT_EPS = param_.getValue("rt_window");
      conf.ENABLE_I = param_.getValue("enable_i").toBool();

      std::string path = param_.getValue("coef_path");
      std::string file = param_.getValue("coef_file");
      std::string coef_path = path + "/" + file;


      PeptideIdentification result = solve_one_prec_mass(spectrum, mass, param_.getValue("rho"), input.getRT(), conf, lastChars, lastCharTags, coef_path);
      
      /** 
        In the following we produce the output for @em input. The individual computations are tagged with a
        tag of the form (n).

        Output Format (on one line):
          OUT;<True Sequence>(1);<Number of results>(1.5);<Best Sequence>(2);<Best Score>(3);...
          ...<Best Predicted RT>(4);<True Score>(5);<True Predicted RT>(6);<True Position>(7);...
          ...<Best Sequence for measured RT>(8);<Best Score for measured RT>(9);...
          ...<Position Best Score for measured RT> (10);ENDOUT;
      */
      LOG("OUTPUT FORMAT: TrueSeq; NrSolutions; BestSeq; BestSeqScore; BestSeqPredRT; TrueSeqScore; TrueSeqPredRT; TrueSeqPos; BestSeqGivenRT; BestSeqGivenRTScore; BestSeqGivenRTPos\n");
      
      // If an annotation sequence is given we search this sequence and output maxscore and score and position of the sequence
      int true_rt = round(input.getRT());
      std::string annotationSequence = param_.getValue("annotation_sequence");
      AASequence annotation;
      std::string not_avail = "-";
      if (annotationSequence.size() > 0)
      {
        // Replace all I with L because they have the same mass and our algorithm only generates L's
        replace(annotationSequence.begin(), annotationSequence.end(), 'I', 'L');
        annotation = AASequence::fromString(annotationSequence);
        std::cout << annotationSequence << ";"; // (1)
      }
      else
        std::cout << not_avail << ";"; // (1)
      
      // We got results
      if (result.getHits().size() > 0)
      {
        std::cout << result.getHits().size() << ";" << result.getHits()[0].getSequence() << ";" << result.getHits()[0].getScore() << ";"; // (1.5) (2) (3)
        PeptideHitRT* phrt = static_cast<PeptideHitRT*>(&result.getHits()[0]);
        if (phrt) 
          std::cout << (int)phrt->getMetaValue("predicted_RT") << ";"; // (4)
        else
          std::cout << not_avail << ";"; // (4)

        // We know the true sequence (annotation)
        if (annotationSequence.size() > 0) 
        {
          // Replace all I with L because they have the same mass and our algorithm only generates L's
          if (!conf.ENABLE_I)
            replace(annotationSequence.begin(), annotationSequence.end(), 'I', 'L');
          int solutionPos = -1;
          for (int i = 0; i < result.getHits().size(); i++) 
          {
            if (result.getHits()[i].getSequence() == annotation) {
              solutionPos = i;
              break;
            }
          }

          if (solutionPos >= 0)
          {
            PeptideHit ph = result.getHits()[solutionPos];
            std::cout << ph.getScore() << ";"; // (5)
          } 
          else
            std::cout << not_avail << ";"; // (5)
          
          if (RetentionTime(coef_path, 0, 0, conf).isNeiModel())
            std::cout << get_seq_rt_nei("-" + annotationSequence + "-") << ";";	// (6)				
          else
            std::cout << get_seq_rt(annotationSequence, 0, 0) << ";";	// (6)				

          if (solutionPos >= 0) 
          {
            while (solutionPos > 0 && result.getHits()[solutionPos-1].getScore() == result.getHits()[solutionPos].getScore())
            {
              solutionPos--;
            }
            std::cout << solutionPos << ";"; // Position in list // (7)
          } 
          else 
            std::cout << not_avail << ";"; // (7)
        }
        // We do not have the true sequence
        else
          std::cout << not_avail << ";"<< not_avail << ";" << not_avail << ";"; // (5) (6) (7)
      } 
      // We do not have any results for this input
      else
        std::cout << 0 << ";" << not_avail << ";"<< not_avail << ";" << not_avail << ";" << not_avail << ";"<< not_avail << ";" << not_avail << ";"; // (1.5) (2) (3) (4) (5) (6) (7)

      // Look if there is a peptide hit that matches the measured RT
      bool found = false;
      for (int i = 0; i < result.getHits().size() && !found; i++) 
      {
        PeptideHitRT* phrt = static_cast<PeptideHitRT*>(&result.getHits()[i]);
        if (phrt && (int)phrt->getMetaValue("predicted_RT") ==  true_rt) 
        {
          std::cout << result.getHits()[i].getSequence() << ";" << result.getHits()[i].getScore() << ";" << i << ";"; // (8) (9) (10)
          found = true;
        }
      }
      if (!found)
        std::cout << not_avail << ";" << not_avail << ";" << not_avail << ";"; // (8) (9) (10)
      std::cout << "ENDOUT" << std::endl;
      
      return result;
    }
  };
}

#endif
