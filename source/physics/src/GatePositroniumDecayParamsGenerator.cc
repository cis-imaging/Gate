/** ----------------------
  Copyright (C): OpenGATE Collaboration
  This software is distributed under the terms
  of the GNU Lesser General  Public Licence (LGPL)
  See LICENSE.md for further details
  ----------------------*/

#include <string>
#include <utility>
#include <vector>

#include "GateMessageManager.hh"
#include "GatePositroniumHelper.hh"
#include "GatePositroniumDecayParamsGenerator.hh"

void GatePositroniumDecayParamsGenerator::SetElectronCaptureProbabilities(const std::vector<float>& electronCaptureProb)
{
  fElectronCaptureProbabilities = electronCaptureProb;
}

void GatePositroniumDecayParamsGenerator::SetPromptGammaProbabilities(const std::vector<float>& promptGammaProb)
{
  fPromptGammaProbabilities = promptGammaProb;
}
void GatePositroniumDecayParamsGenerator::SetDecayKinds(const std::vector<PositroniumDecayKind>& decayKinds)
{
  fDecayKinds = decayKinds;
}
void GatePositroniumDecayParamsGenerator::SetPositronInteractions(const std::vector<PositronElectronInteraction>& positronInteractions)
{
  fPositronInteractions = positronInteractions;
}
void GatePositroniumDecayParamsGenerator::SetPositroniumLifetimes(const std::vector<float>& positroniumLifetimes)
{
  fPositroniumLifetimes = positroniumLifetimes;
}
void GatePositroniumDecayParamsGenerator::SetPromptGammaEnergies(const std::vector<float>& energies)
{
  fPromptGammaEnergies = energies;
}
void GatePositroniumDecayParamsGenerator::SetPositroniumFraction(const std::vector<float>& positroniumFractions)
{
  fPositroniumFractions = positroniumFractions;
}

void GatePositroniumDecayParamsGenerator::SetMeanPositronRange(const std::vector<float>& meanPositironRange)
{
  fMeanPositronRange = meanPositironRange;
}


void GatePositroniumDecayParamsGenerator::RaiseErrorIfModelParamsNotSet() const
{
  const std::vector<std::pair<std::string, bool>> paramsSet = {
    {"Decay kinds", fDecayKinds.has_value()},
    {"Positron interactions", fPositronInteractions.has_value()},
    {"Positronium decay fractions", fPositroniumFractions.has_value()},
    {"Positronium lifetimes", fPositroniumLifetimes.has_value()},
    {"Prompt gamma probabilities", fPromptGammaProbabilities.has_value()},
    {"Prompt gamma energies", fPromptGammaEnergies.has_value()}
  };
  for (const auto& [name, isSet] : paramsSet) {
    if (!isSet) {
      GateError("GatePositroniumDecayParamsGenerator::RaiseErrorIfModelParamsNotSet: "
                << name << " are not set");
    }
  }
}

void GatePositroniumDecayParamsGenerator::RaiseErrorIfModelParamsSizesDiffer() const
{
  // The number of channels is defined by setDecayKinds; every other vector must match it.
  const std::size_t nChannels = fDecayKinds->size();
  if (nChannels == 0) {
    GateError("GatePositroniumDecayParamsGenerator::RaiseErrorIfModelParamsSizesDiffer: "
              "Decay kinds are empty");
  }

  std::vector<std::pair<std::string, std::size_t>> paramsSize = {
    {"Positron interactions", fPositronInteractions->size()},
    {"Positronium decay fractions", fPositroniumFractions->size()},
    {"Positronium lifetimes", fPositroniumLifetimes->size()},
    {"Prompt gamma probabilities", fPromptGammaProbabilities->size()},
    {"Prompt gamma energies", fPromptGammaEnergies->size()}
  };
  if (fElectronCaptureProbabilities.has_value()) {
    paramsSize.emplace_back("Electron capture probabilities", fElectronCaptureProbabilities->size());
  }
  if (fMeanPositronRange.has_value()) {
    paramsSize.emplace_back("Mean positron ranges", fMeanPositronRange->size());
  }

  for (const auto& [name, size] : paramsSize) {
    if (size != nChannels) {
      GateError("GatePositroniumDecayParamsGenerator::RaiseErrorIfModelParamsSizesDiffer: "
                << name << " have " << size << " elements, but decay kinds define "
                << nChannels << " channels");
    }
  }
}

void GatePositroniumDecayParamsGenerator::RaiseErrorIfChannelsUnphysical() const
{
  // Para-positronium annihilates into an even number of gammas, so a kParaPs channel cannot
  // be k3Gamma. kOrthoPs + k2Gamma (pick-off, quenching) and kDirect + k3Gamma are allowed.
  for (std::size_t i = 0; i < fDecayKinds->size(); ++i) {
    if ((*fPositronInteractions)[i] == PositronElectronInteraction::kParaPs &&
        (*fDecayKinds)[i] == PositroniumDecayKind::k3Gamma) {
      GateError("GatePositroniumDecayParamsGenerator::RaiseErrorIfChannelsUnphysical: "
                "channel " << i << " is kParaPs with k3Gamma; para-positronium decays into 2 gammas");
    }
  }
}

void GatePositroniumDecayParamsGenerator::RaiseErrorIfFractionsInvalid() const
{
  // Fractions are relative weights normalized by GatePositroniumHelper::NormalizeFractions.
  // A negative weight gives probabilities outside [0, 1]; a zero sum cannot be normalized and
  // makes the model always pick the last channel.
  double sum = 0.0;
  for (std::size_t i = 0; i < fPositroniumFractions->size(); ++i) {
    const float fraction = (*fPositroniumFractions)[i];
    if (fraction < 0.0f) {
      GateError("GatePositroniumDecayParamsGenerator::RaiseErrorIfFractionsInvalid: "
                "fraction of channel " << i << " is negative (" << fraction << ")");
    }
    sum += fraction;
  }
  if (sum <= 0.0) {
    GateError("GatePositroniumDecayParamsGenerator::RaiseErrorIfFractionsInvalid: "
              "sum of positronium decay fractions must be positive");
  }
}

PositroniumDecayModelParams GatePositroniumDecayParamsGenerator::generatePositroniumDecayParams() const
{
  RaiseErrorIfModelParamsNotSet();
  RaiseErrorIfModelParamsSizesDiffer();
  RaiseErrorIfChannelsUnphysical();
  RaiseErrorIfFractionsInvalid();

  // Every user-defined component is exactly one channel of the model (setDecayKinds and
  // setPositronInteractions are both required). Nothing is expanded here, so all vectors
  // have nChannels elements and are indexed by the same decayIndex.
  const std::size_t nChannels = fDecayKinds->size();
  GatePositroniumHelper positronHelper;

  PositroniumDecayModelParams params;
  params.fDecayKind = *fDecayKinds;
  params.fPositronInteractions = *fPositronInteractions;
  params.fFractions = positronHelper.NormalizeFractions(*fPositroniumFractions);
  params.fLifetimes = *fPositroniumLifetimes;
  params.fPromptGammaProbabilities = *fPromptGammaProbabilities;
  params.fPromptGammaEnergy = *fPromptGammaEnergies;
  params.fElectronCaptureProbabilities =
    fElectronCaptureProbabilities.value_or(std::vector<float>(nChannels, 0.0f));
  params.fMeanPositronRange = fMeanPositronRange.value_or(std::vector<float>(nChannels, 0.0f));
  params.fMeanPositronRangeEnabled.assign(nChannels, fMeanPositronRange.has_value());

  validatePositroniumDecayParams(params);
  return params;
}

void GatePositroniumDecayParamsGenerator::validatePositroniumDecayParams(const PositroniumDecayModelParams& params) const
{
  // Invariant of the final structure: GatePositroniumDecayModel indexes each of these
  // vectors by the channel index, so all of them must have as many elements as fDecayKind.
  const std::size_t nChannels = params.fDecayKind.size();
  const bool sizeMismatch = (params.fFractions.size() != nChannels) ||
                            (params.fLifetimes.size() != nChannels) ||
                            (params.fPositronInteractions.size() != nChannels) ||
                            (params.fPromptGammaProbabilities.size() != nChannels) ||
                            (params.fPromptGammaEnergy.size() != nChannels) ||
                            (params.fElectronCaptureProbabilities.size() != nChannels) ||
                            (params.fMeanPositronRangeEnabled.size() != nChannels) ||
                            (params.fMeanPositronRange.size() != nChannels);
  if (sizeMismatch) {
    GateError("GatePositroniumDecayParamsGenerator::validatePositroniumDecayParams: "
              "not every per-channel parameter vector has " << nChannels << " elements");
  }

  constexpr double kProb_of_no_particle_limit = 0.9;
  for (std::size_t i = 0; i < nChannels; ++i) {
    const double electronCaptureProb = params.fElectronCaptureProbabilities[i];
    const double noPromptProb = 1.0 - params.fPromptGammaProbabilities[i];
    if (electronCaptureProb * noPromptProb >= kProb_of_no_particle_limit) {
      GateWarning("The probability of 0 prompt emission times probability of electron capture (no anihillation) is larger than 90%. Please check channel definitions! Otherwise simulations can take a lot of time");
    }
  }
}
