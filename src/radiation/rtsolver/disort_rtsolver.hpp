// Athena++ headers
#include "../../athena.hpp"

class RadiationBand;
class Radiation;

struct disort_state;
struct disort_output;

class DisortRTSolver: public RTSolver {
public:
  RadiationBand *pmy_band;

  DisortRTSolver(RadiationBand *pband, ParameterInput *pin):
    RTSolver(pband, pin) {};
  ~DisortRTSolver();

  void RadiativeTransfer(int n, int k, int j, int il, int iu);

protected:
  disort_state *ds;
  disort_output *ds_out;

  void init_disort(ParameterInput *pin);
  void free_disort();
};