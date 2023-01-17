// Athena++ headers
#include "../../athena.hpp"
#include "../../globals.hpp"
#include "rtsolver.hpp"


class RadiationBand;
class Radiation;

struct disort_state;
struct disort_output;

class DisortRTSolver: public RTSolver {
public:
  RadiationBand *pmy_band;

  DisortRTSolver(RadiationBand *pband, ParameterInput *pin);
  ~DisortRTSolver();

  void RadiativeTransfer(MeshBlock *pmb, int n, int k, int j);

protected:
  disort_state *ds;
  disort_output *ds_out;

  void init_disort(ParameterInput *pin);
  void free_disort();
};