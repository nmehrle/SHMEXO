// C++ headers
#include <sstream>

// Athena++ headers
#include "hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../globals.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../utils/buffer_utils.hpp"
#include "vertical_communication.hpp"

int TAG_TOPPRESSURE = 1111;
int TAG_BOTPRESSURE = 2222;

Hydro::~Hydro() {
  delete psbuf_;
  delete pvc;
}

void Hydro::CheckHydro() {
  MeshBlock *pmb = pmy_block->pmy_mesh->pblock;
  std::stringstream msg;
  int myrank = Globals::my_rank;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        for (int n = 0; n < NMASS; ++n)
          if (w(n,k,j,i) < 0.) {
            msg << "### FATAL ERROR in Hydro::CheckHydro" << std::endl
                << "Density variable is negative at position ("
                << k << "," << j << "," << i << ") in rank " << myrank;
            ATHENA_ERROR(msg);
          }
        if (w(IPR,k,j,i) < 0.) {
          msg << "### FATAL ERROR in Hydro::CheckHydro" << std::endl
              << "Pressure is negative at position ("
              << k << "," << j << "," << i << ") in rank " << myrank;
          ATHENA_ERROR(msg);
        }
        Real RT = w(IPR,k,j,i)/w(IDN,k,j,i);
        Real grav = -hsrc.GetG1(k, j, i);
        if (grav != 0) {
          Real RTmin = 2.*grav*pmb->pcoord->dx1f(i);
          if (RT < RTmin) {
            msg << "### FATAL ERROR in Hydro::CheckHydro" << std::endl
                << "Vertical spacing is less than half scale height at position ("
                << k << "," << j << "," << i << ") in rank " << myrank << std::endl;
            ATHENA_ERROR(msg);
          }
        }
      }

  // make a copy of w, needed for outflow boundary condition
  w1 = w;
  if (myrank == 0)
    std::cout << "Hydro check passed. ";
}

// FIXME: local boundary has not been implemented
// Ordering the meshblocks need to be worked out such that
// the upper boundary executes before the lower boundary
void Hydro::RecvTopPressure(AthenaArray<Real> &psf, AthenaArray<Real> &entropy,
  AthenaArray<Real> &gamma, NeighborBlock nbot,
  int kl, int ku, int jl, int ju)
{
  int ie = pmy_block->ie;
  int ssize = 3*(ju-jl+1)*(ku-kl+1);

  std::stringstream msg;
#ifdef MPI_PARALLEL
  MPI_Status status;
#endif

  if (nbot.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
    int tag = BoundaryBase::CreateBvalsMPITag(pmy_block->lid, TAG_TOPPRESSURE, nbot.bufid);
    MPI_Recv(psbuf_, ssize, MPI_ATHENA_REAL, nbot.snb.rank, tag, MPI_COMM_WORLD, &status);
#endif
  } else {  // local boundary
    // need to wait for the top boundary to finish
    msg << "### FATAL ERROR in Hydro::RecvTopPressure" << std::endl
        << "Local boundary not yet implemented" << std::endl;
    ATHENA_ERROR(msg);
  }
  int p = 0;
  BufferUtility::UnpackData(psbuf_, psf, ie+1, ie+1, jl, ju, kl, ku, p);
  BufferUtility::UnpackData(psbuf_, entropy, ie, ie, jl, ju, kl, ku, p);
  BufferUtility::UnpackData(psbuf_, gamma, ie, ie, jl, ju, kl, ku, p);
}

void Hydro::SendTopPressure(AthenaArray<Real> &psf, AthenaArray<Real> &entropy,
  AthenaArray<Real> &gamma, NeighborBlock ntop,
  int kl, int ku, int jl, int ju)
{
  int is = pmy_block->is, ie = pmy_block->ie;
  int ssize = 0;

  BufferUtility::PackData(psf, psbuf_, is, is, jl, ju, kl, ku, ssize);
  BufferUtility::PackData(entropy, psbuf_, ie, ie, jl, ju, kl, ku, ssize);
  BufferUtility::PackData(gamma, psbuf_, ie, ie, jl, ju, kl, ku, ssize);
  if (ntop.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
    int tag = BoundaryBase::CreateBvalsMPITag(ntop.snb.lid, TAG_TOPPRESSURE, ntop.targetid);
    MPI_Isend(psbuf_, ssize, MPI_ATHENA_REAL, ntop.snb.rank, tag, MPI_COMM_WORLD,
      &req_send_top_pressure_);
#endif
  } else {  // local boundary
    MeshBlock *pbl = pmy_block->pmy_mesh->FindMeshBlock(ntop.snb.gid);
    std::memcpy(pbl->phydro->psbuf_, psbuf_, ssize*sizeof(Real));
  }
}

void Hydro::WaitTopPressure()
{
#ifdef MPI_PARALLEL
  MPI_Status status;
  MPI_Wait(&req_send_top_pressure_, &status);
#endif
}

// Ordering the meshblocks should ensure that
// the lower boundary executes before the upper boundary
void Hydro::RecvBotPressure(AthenaArray<Real> &psf, AthenaArray<Real> &entropy,
  AthenaArray<Real> &gamma, NeighborBlock nbot,
  int kl, int ku, int jl, int ju)
{
  int is = pmy_block->is;
  int ssize = 3*(ju-jl+1)*(ku-kl+1);

  std::stringstream msg;
#ifdef MPI_PARALLEL
  MPI_Status status;
#endif

  if (nbot.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
    int tag = BoundaryBase::CreateBvalsMPITag(pmy_block->lid, TAG_BOTPRESSURE, nbot.bufid);
    MPI_Recv(psbuf_, ssize, MPI_ATHENA_REAL, nbot.snb.rank, tag, MPI_COMM_WORLD, &status);
#endif
  } // local boundary will automatically work because it is treated in SendBotPressure

  int p = 0;
  BufferUtility::UnpackData(psbuf_, psf, is, is, jl, ju, kl, ku, p);
  BufferUtility::UnpackData(psbuf_, entropy, is, is, jl, ju, kl, ku, p);
  BufferUtility::UnpackData(psbuf_, gamma, is, is, jl, ju, kl, ku, p);
}

void Hydro::SendBotPressure(AthenaArray<Real> &psf, AthenaArray<Real> &entropy,
  AthenaArray<Real> &gamma, NeighborBlock ntop,
  int kl, int ku, int jl, int ju)
{
  int is = pmy_block->is, ie = pmy_block->ie;
  int ssize = 0;

  BufferUtility::PackData(psf, psbuf_, ie+1, ie+1, jl, ju, kl, ku, ssize);
  BufferUtility::PackData(entropy, psbuf_, is, is, jl, ju, kl, ku, ssize);
  BufferUtility::PackData(gamma, psbuf_, is, is, jl, ju, kl, ku, ssize);
  if (ntop.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
    int tag = BoundaryBase::CreateBvalsMPITag(ntop.snb.lid, TAG_BOTPRESSURE, ntop.targetid);
    MPI_Isend(psbuf_, ssize, MPI_ATHENA_REAL, ntop.snb.rank, tag, MPI_COMM_WORLD,
      &req_send_bot_pressure_);
#endif
  } else {  // local boundary
    MeshBlock *pbl = pmy_block->pmy_mesh->FindMeshBlock(ntop.snb.gid);
    std::memcpy(pbl->phydro->psbuf_, psbuf_, ssize*sizeof(Real));
  }
}

void Hydro::WaitBotPressure()
{
#ifdef MPI_PARALLEL
  MPI_Status status;
  MPI_Wait(&req_send_bot_pressure_, &status);
#endif
}
