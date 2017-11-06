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

#include <OpenMS/METADATA/PeptideHit.h>

#ifndef OPENMS_ANALYSIS_DENOVO_SYMDIFF_PEPTIDEHITRT_H
#define OPENMS_ANALYSIS_DENOVO_SYMDIFF_PEPTIDEHITRT_H

/**
  @brief Class that extends the class @em PeptideHit with a field to store the retention time.
*/
class PeptideHitRT : public PeptideHit {
protected:	
  const int pred_rt_;
public:
  PeptideHitRT(double score, UInt rank, Int charge, const AASequence& sequence, int predicted_rt) : 
    PeptideHit(score, rank, charge, sequence), pred_rt_(predicted_rt)
    {
      PeptideHit::setMetaValue("predicted_RT", pred_rt_);
    }

  int getPredRT() const { return pred_rt_; }

  // assignment operator
  PeptideHitRT& operator=(const PeptideHitRT& source);

  // Equality operator
  bool operator==(const PeptideHitRT& rhs) const;

  // Inequality operator
  bool operator!=(const PeptideHitRT& rhs) const;

};

bool PeptideHitRT::operator==(const PeptideHitRT& rhs) const
{
  return PeptideHit::operator==(rhs) && pred_rt_ == rhs.pred_rt_;	
}

bool PeptideHitRT::operator!=(const PeptideHitRT& rhs) const
{
  return !operator==(rhs);
}

#endif
