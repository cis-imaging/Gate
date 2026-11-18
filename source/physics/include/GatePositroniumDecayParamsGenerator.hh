/** ----------------------
  Copyright (C): OpenGATE Collaboration
  This software is distributed under the terms
  of the GNU Lesser General  Public Licence (LGPL)
  See LICENSE.md for further details
  ----------------------*/

/** Authors: Wojciech Krzemień, Mateusz Bała and Kamil Dulski
 *  Emails: wojciech.krzemien@ncbj.gov.pl, mateusz.bala@ncbj.gov.pl and kamil.dulski@gmail.com
 *  Organization: National Centre For Nuclear Research (NCBJ, https://ncbj.gov.pl), Poland
 *  Developed within the IMPET project: https://pet.ncbj.gov.pl/
 *  About class: Builder that assembles and validates a PositroniumDecayModelParams from individually set parameters. Every user-defined component is exactly one decay channel of the model.
 **/

#ifndef GatePositroniumDecayParamsGenerator_hh
#define GatePositroniumDecayParamsGenerator_hh

#include <optional>
#include <vector>

#include "GatePositroniumDecayModelParams.hh"

/*! class GatePositroniumDecayParamsGenerator
 *  brief Builds PositroniumDecayModelParams for GatePositroniumDecayModel.
 *
 *  Decay kinds, positron interactions, fractions, lifetimes, prompt-gamma probabilities and
 *  prompt-gamma energies are required. Electron capture probabilities (default 0) and mean
 *  positron ranges (default: range disabled) are optional. All vectors must have one entry
 *  per decay kind; a missing required parameter or a length mismatch is reported with GateError.
 *  A kParaPs channel must be k2Gamma. Fractions must be non-negative with a positive sum.
 */


class GatePositroniumDecayParamsGenerator
{
public:

  GatePositroniumDecayParamsGenerator()=default;
  virtual ~GatePositroniumDecayParamsGenerator()=default;

  void SetElectronCaptureProbabilities(const std::vector<float>& electronCaptureProb);
  void SetPromptGammaProbabilities(const std::vector<float>& promptGammaProb);
  void SetPromptGammaEnergies(const std::vector<float>& energies);
  void SetPositroniumLifetimes(const std::vector<float>& fPositroniumLifetimes);
  void SetDecayKinds(const std::vector<PositroniumDecayKind>& decayKinds);
  void SetPositronInteractions(const std::vector<PositronElectronInteraction>& positronInteractions);
  void SetPositroniumFraction(const std::vector<float>& positroniumFractions);
  void SetMeanPositronRange(const std::vector<float>& meanPositronRange);

  PositroniumDecayModelParams generatePositroniumDecayParams() const;
  void validatePositroniumDecayParams(const PositroniumDecayModelParams& params) const;


private:
  void RaiseErrorIfModelParamsNotSet() const;
  void RaiseErrorIfModelParamsSizesDiffer() const;
  void RaiseErrorIfChannelsUnphysical() const;
  void RaiseErrorIfFractionsInvalid() const;

  std::optional<std::vector<float>> fPositroniumFractions;
  std::optional<std::vector<float>> fPositroniumLifetimes;
  std::optional<std::vector<float>> fPromptGammaProbabilities;
  std::optional<std::vector<float>> fPromptGammaEnergies;
  std::optional<std::vector<float>> fElectronCaptureProbabilities;
  std::optional<std::vector<PositroniumDecayKind>> fDecayKinds;
  std::optional<std::vector<PositronElectronInteraction>> fPositronInteractions;
  std::optional<std::vector<float>> fMeanPositronRange;
};
#endif
