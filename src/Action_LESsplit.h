#ifndef INC_ACTION_LESSPLIT_H
#define INC_ACTION_LESSPLIT_H
#include "Action.h"
#include "EnsembleOut_Multi.h"
#include "Trajout_Single.h"
/// Split LES frame/top into normal frames/tops.
class Action_LESsplit : public Action {
  public:
    Action_LESsplit() : lesAverage_(false), lesSplit_(false), lesParm_(0) {}
    ~Action_LESsplit();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_LESsplit(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    bool lesAverage_;
    bool lesSplit_;
    typedef std::vector<AtomMask> MaskArray;
    MaskArray lesMasks_;
    EnsembleOut_Multi lesTraj_;
    Trajout_Single avgTraj_;
    std::string trajfilename_;
    std::string avgfilename_;
    ArgList trajArgs_;
    FrameArray lesFrames_;
    FramePtrArray lesPtrs_;
    Frame avgFrame_;
    Topology* lesParm_;
};
#endif
