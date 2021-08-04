#ifndef VERTICAL_COMMUNICATION_HPP_
#define VERTICAL_COMMUNICATION_HPP_

// C/C++ headers
#include <cstring>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../bvals/bvals_interfaces.hpp"
#include "../mesh/mesh.hpp"
#include "hydro.hpp"

// MPI headers
#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

class VerticalCommunication {
public:
// data
  bool has_top_neighbor, has_bot_neighbor;
  NeighborBlock tblock, bblock;

// functions
  VerticalCommunication(Hydro *phydro);
  ~VerticalCommunication();
  void WaitSendTopImplicit(int kl, int ku, int jl, int ju);
  void WaitSendBotImplicit(int kl, int ku, int jl, int ju);
  int CreateMPITag(int lid, int bufid, int phy);
  void FindNeighbors();

  template<typename T1, typename T2>
  void RecvBotImplicitData(T1 &a, T2 &b, int k, int j, NeighborBlock nbot) {
    int s1 = a.size(), s2 = b.size();
    int phy = k << 6 | j << 1 | 1;
#ifdef MPI_PARALLEL
    MPI_Status status;
#endif

    if (nbot.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
      int tag = CreateMPITag(pmy_hydro_->pmy_block->lid, nbot.bufid, phy);
      MPI_Recv(buffer_, s1+s2, MPI_ATHENA_REAL, nbot.snb.rank, tag, MPI_COMM_WORLD, &status);
#endif
    } // local boundary

    memcpy(a.data(), buffer_, s1*sizeof(Real));
    memcpy(b.data(), buffer_ + s1, s2*sizeof(Real));
  }

  template<typename T1, typename T2>
  void SendTopImplicitData(T1 &a, T2 &b, int k, int j, NeighborBlock ntop) {
    int s1 = a.size(), s2 = b.size();
    int phy = k << 6 | j << 1 | 1;

    memcpy(buffer_, a.data(), s1*sizeof(Real));
    memcpy(buffer_ + s1, b.data(), s2*sizeof(Real));

    if (ntop.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
      int tag = CreateMPITag(ntop.snb.lid, ntop.targetid, phy);
      MPI_Isend(buffer_, s1+s2, MPI_ATHENA_REAL, ntop.snb.rank, tag, MPI_COMM_WORLD,
        &req_send_top_data_[k][j]);
#endif
    } // local boundary
  }

  template<typename T>
  void RecvTopImplicitData(T &a, int k, int j, NeighborBlock ntop) {
    int s1 = a.size();
    int phy = k << 6 | j << 1 | 0;
#ifdef MPI_PARALLEL
    MPI_Status status;
#endif

    if (ntop.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
      int tag = CreateMPITag(pmy_hydro_->pmy_block->lid, ntop.bufid, phy);
      MPI_Recv(buffer_, s1, MPI_ATHENA_REAL, ntop.snb.rank, tag, MPI_COMM_WORLD, &status);
#endif
    } // local boundary

    memcpy(a.data(), buffer_, s1*sizeof(Real));
  }

  template<typename T>
  void SendBotImplicitData(T &a, int k, int j, NeighborBlock nbot) {
    int s1 = a.size();
    int phy = k << 6 | j << 1 | 0;

    memcpy(buffer_, a.data(), s1*sizeof(Real));

    if (nbot.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
      int tag = CreateMPITag(nbot.snb.lid, nbot.targetid, phy);
      MPI_Isend(buffer_, s1, MPI_ATHENA_REAL, nbot.snb.rank, tag, MPI_COMM_WORLD,
        &req_send_bot_data_[k][j]);
#endif
    } // local boundary
  }

  template<typename T1, typename T2>
  void SaveImplicitData(std::vector<T1> &a, std::vector<T2> &b,
    int k, int j, int il, int iu) {
    for (int i = il; i <= iu; ++i) {
      int s1 = a[i].size(), s2 = b[i].size();
      memcpy(implicit_data_[k][j][i], a[i].data(), s1*sizeof(Real));
      memcpy(implicit_data_[k][j][i] + s1, b[i].data(), s2*sizeof(Real));
    }
  }

  template<typename T1, typename T2>
  void LoadImplicitData(std::vector<T1> &a, std::vector<T2> &b,
    int k, int j, int il, int iu) {
    for (int i = il; i <= iu; ++i) {
      int s1 = a[i].size(), s2 = b[i].size();
      memcpy(a[i].data(), implicit_data_[k][j][i], s1*sizeof(Real));
      memcpy(b[i].data(), implicit_data_[k][j][i] + s1, s2*sizeof(Real));
    }
  }

private:
  Hydro *pmy_hydro_;
  Real *buffer_;
  Real ****implicit_data_; // archive of implicit tri-diagonal matrix

#ifdef MPI_PARALLEL
  MPI_Request **req_send_bot_data_;
  MPI_Request **req_send_top_data_;
#endif
};

#endif
