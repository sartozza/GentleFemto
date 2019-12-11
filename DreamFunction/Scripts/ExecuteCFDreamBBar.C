#include "GetCorrelationsNano.C"
#include "GetCorrelationsNanoBBar.C"
#include "GetCorrelationsBbarB.C"
#include "GetCorrelationsNanoBB.C"
#include "METoSEReweighting.C"
#include "TSystem.h"
int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
   GetCorrelationsNanoBBar(filename, prefix, addon, 0.18, 0.28);
  //  GetCorrelationsNanoBB(filename, prefix, addon);

  return 1;
}
