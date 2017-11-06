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

#ifndef OPENMS_ANALYSIS_DENOVO_SYMDIFF_DATATYPES_H
#define OPENMS_ANALYSIS_DENOVO_SYMDIFF_DATATYPES_H

#include <stdio.h>
#include <vector>
#include <string>
#include <map>

#define mp make_pair
// Lightweight container for peptides; use as follows ((<score>, <retention time>), <amino acid sequence>)
#define peptide std::pair<std::pair<int, int>, std::vector<std::string> >

// @em edge is an edge of our matrix spectrum graph
typedef struct 
{
  int from_x;
  int from_y;
  int to_x;
  int to_y;
  std::string label;
  std::string label_x;
  std::string label_y;
} edge;

// @em peak is a simple peak of a mass spectrum
typedef struct 
{
  double time;
  double mz;
  double h;
} peak;

// Sorting function to sort peaks by their time and mass
struct comp_peak 
{
  bool operator() (const peak& lhs, const peak& rhs) const 
  {
    if (lhs.time != rhs.time) 
      return lhs.time < rhs.time;
    return lhs.mz < rhs.mz;
  }
};

#endif
