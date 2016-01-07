#ifndef INC_ACTION_GIBBSEOFHYDRATION_H
#define INC_ACTION_GIBBSEOFHYDRATION_H
#include "Action.h"
#include "DataSet_GridFlt.h"
#include "GridAction.h"
/** \author Mark J. Williamson
  * \author C++ adaptation by DRR
  *  For more on the theory, please see eq. 1 in http://dx.doi.org/10.1021/ci100462t 
  */
class Action_GridFreeEnergy : public Action, private GridAction {
  public:
    Action_GridFreeEnergy();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_GridFreeEnergy(); }
    void Help() const;
  private:
    // Action members
    Action::RetType Init(ArgList&, ActionInit&, int);
#   ifdef MPI
    int ParallelActionInit(Parallel::Comm const& c) { return ParallelGridInit(c, grid_); }
#   endif
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();

    /// maximum expected voxel occupancy count
    // TODO Work out a smart way of calculating this since it will be a 
    //      function of the trajectory.Maybe an upper limit of this is:
    //          numberOfVoxels * numberOfFrames?
    int maxVoxelOccupancyCount_;
    // Temperature to calculate gfe at
    double tempInKevin_;
    /// Atom mask
    AtomMask mask_;
    DataSet_GridFlt* grid_;
};
#endif
