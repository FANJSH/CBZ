#ifndef PROMPTKRATIO
#define PROMPTKRATIO

#include "Numeric.h"
#include "LibData.h"
#include "DelayedNeutronData.h"
#include "SNR_system.h"
#include "SensitivityData.h"

class PromptKratio{
 protected:
 public:
  PromptKratio();
  void DelayedNeutronDataExtractionFromChi
    (XSLibrary& xslib, DelayedNeutronData& dnd, int matid);
  void DelayedNeutronDataAdditionToChi
    (XSLibrary& xslib, DelayedNeutronData& dnd, int matid, real factor);
  SensitivityData CalSensitivity
    //    (SNRSystem* fwd, SNRSystem* adj, XSLibrary &xslib, DelayedNeutronData& dnd,
    (GeneralSystem* fwd, GeneralSystem* adj, XSLibrary &xslib, DelayedNeutronData& dnd,
     real keff, int nucnum, int *nucid, real factor);
};

#endif
