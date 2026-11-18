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
 *  About class: Utility class computing per-component decay fractions from lifetimes and positron-electron interaction types (ortho-Ps, direct annihilation, para-Ps), using the hyperfine coefficient and mean ortho-Ps lifetime; also normalizes fraction vectors.
 **/

#ifndef GatePositroniumHelper_h
#define GatePositroniumHelper_h 1

#include "GatePositroniumDecayModelParams.hh"

class GatePositroniumHelper
{
public:
  GatePositroniumHelper() {};
  virtual ~GatePositroniumHelper() {};

  /// Expands each non-pPs component into a 2-gamma and a 3-gamma channel, with intensities
  /// derived from the lifetimes, and appends a pPs channel.
  /// NOTE: currently not used anywhere. GatePositroniumDecayParamsGenerator requires
  /// setDecayKinds and maps every user component to exactly one channel. Kept for possible
  /// future use; CalcFractionsFromLifetime and CalcPPsFractionFromOPs are only called from here.
  /// The function expands neither fMeanPositronRangeEnabled nor, when it is empty,
  /// fMeanPositronRange - fix that before using it again.
  struct PositroniumDecayModelParams CalculateFractionsFromLifetimes(struct PositroniumDecayModelParams params);
  float CalcPPsFractionFromOPs(std::vector<float> fractions, std::vector<PositronElectronInteraction> decays);
  std::pair<float, float> CalcFractionsFromLifetime(float intensity, float lifetime, PositronElectronInteraction inter);
  float CalcFractionFromOPsLifetime(float intensity, float lifetime, PositroniumDecayKind decay);
  float CalcFractionFromDirectLifetime(float intensity, PositroniumDecayKind decay);
  std::vector<float> NormalizeFractions (std::vector<float> fractions);
};

#endif
