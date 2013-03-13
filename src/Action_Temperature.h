#ifndef INC_ACTION_TEMPERATURE_H
#define INC_ACTION_TEMPERATURE_H
#include "Action.h"
/// Calculate the temperature of parts of a system.
class Action_Temperature : public Action {
  public:
    Action_Temperature();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Temperature(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataSet* Tdata_;
    AtomMask Mask_;
    int ntc_;
    int degrees_of_freedom_;
};
#endif
