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

#ifndef OPENMS_ANALYSIS_DENOVO_SYMDIFF_CONFIG_H
#define OPENMS_ANALYSIS_DENOVO_SYMDIFF_CONFIG_H

#include <OpenMS/DATASTRUCTURES/DataValue.h>

/**
  @em config is a struct that holds the configuration parameters for the preprocessing of the the spectrum
  and for our algorithm.
*/ 
typedef struct {
  double EPS; // Accuracy used for equality of peaks (in Daltons)
  double MERGING_EPS; // Accuracy used for merging peaks (in Daltons)
  double NOISE_H_CUTOFF; // Peaks with lower intensity are ignored
  double MISSING_PEAK_PUNISHMENT; // Negative score for a missing b-/y-ion peak pair
  double HALF_MISSING_PEAK_PUNISHMENT; // Punishment if only one peak of a b-/y-ion peak pair is missing
  double MISSING_PEAK_THRESHOLD_SCORE; // Maximum intensity a peak can have to be counted as missing
  double MISSING_ION_PEAK_PUNISHMENT; // Punishment if a ion other than b-/y-ion is missing
  double START_SCORE; // Score for an empty path pair
  bool COUNT_PEAKS; // Count number of peaks instead of summing up their intensities
  OpenMS::DoubleList b_ion_offsets; // Offsets for additional ion types to the b-ion
  OpenMS::DoubleList y_ion_offsets; // Offsets for additional ion types to the y-ion

  // Parameters for RT prediction
  bool ENABLE_I; // Consider 'I' and 'L' amino acids separately; their masses are equal but their RT coefficients not
  int RT_EPS; // Retention time accuracy
} config;

#endif
