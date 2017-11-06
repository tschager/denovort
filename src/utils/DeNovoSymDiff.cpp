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

#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/DENOVO/SYMDIFF/DeNovoSymDiffAlgorithm.h>
#include <sstream>

using namespace OpenMS;


class DeNovoSymDiff
  : public TOPPBase
{
 public:

  DeNovoSymDiff()
    : TOPPBase("DeNovoSymDiff","Finds a amino acid sequence for a given spectrum  and a given retention time.", false)
  {
  }

 protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in","<file>","","input file ");
    setValidFormats_("in",ListUtils::create<String>("dta2d"));

    registerOutputFile_("out","<file>","","output file");
    setValidFormats_("out",ListUtils::create<String>("idXML"));

    registerDoubleOption_("parent_mass","<value>",0,"mass of the peptide (in Daltons)", false);

    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String &) const
  {
    return DeNovoSymDiffAlgorithm().getDefaults();
  }

  ExitCodes main_(int , const char**)
  {

    String in = getStringOption_("in");
    String out = getStringOption_("out");

    Param params = getParam_().copy("algorithm:", true);

    MSExperiment<Peak1D> exp;

    DTA2DFile f;
    f.setLogType(log_type_);
    f.load(in,exp);

    DeNovoSymDiffAlgorithm sdd;
    sdd.setParameters(params);

    std::vector<PeptideIdentification> pep_ids;

    for(Size i=0; i<exp.size(); i++)
    {
        pep_ids.push_back(sdd.searchSequences(exp[i], getDoubleOption_("parent_mass")));
    }


    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    DateTime now = DateTime::now();
    String date_string = now.get();
    String identifier("CompNovo_" + date_string);

    for (std::vector<PeptideIdentification>::iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
    {
      it->setIdentifier(identifier);
    }

    std::vector<ProteinIdentification> prot_ids;
    ProteinIdentification prot_id;
    prot_id.setIdentifier(identifier);
    prot_id.setDateTime(now);
    prot_id.setSearchEngineVersion("0.1");
    prot_id.setSearchEngine("DeNovoSymDiff");
    prot_ids.push_back(prot_id);

    IdXMLFile().store(out, prot_ids, pep_ids);

    return EXECUTION_OK;
  }

};

int main( int argc, const char** argv )
{
  DeNovoSymDiff tool;
  return tool.main(argc,argv);
}

