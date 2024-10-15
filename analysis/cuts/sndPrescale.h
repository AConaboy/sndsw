#pragma once

#include "sndBaseCut.h"

#include <random>

namespace snd{
  namespace analysis_cuts{
    class preScale : public snd::analysis_cuts::baseCut{
    private:

      std::random_device r;
      std::mt19937 gen;
      std::uniform_real_distribution<double> dis;

      double fraction_;
      
    public:
      preScale(double fraction);
      bool passCut();
    };
  }
}
