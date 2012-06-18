#ifndef INC_ACTION_HBOND_H
#define INC_ACTION_HBOND_H
#include <vector>
#include <map>
#include <set>
#include "Action.h"
// Class: Action_Hbond
/// Action to calculate the Hbonds present in each frame.
class Action_Hbond : public Action {
  public:
    Action_Hbond();

    void print();
  private:
    int init();
    int setup();
    int action();

    struct HbondType {
      int A;        ///< Acceptor atom#
      int H;        ///< Hydrogen atom#
      int D;        ///< Donor atom#
      int Frames;   ///< # frames this hbond has been present
      double dist;  ///< Used to calc avg distance of this hbond
      double angle; ///< Used to calc avg angle of this hbond
    };

    int Nframes_;
    std::string avgout_;
    std::string solvout_;
    std::string bridgeout_;
    typedef std::map<int,HbondType> HBmapType;
    HBmapType HbondMap_;   ///< Track all solute-solute hbonds found
    HBmapType SolventMap_; ///< Track all solute-solvent hbonds found
    typedef std::map< std::set<int>, int > BridgeType;
    BridgeType BridgeMap_; ///< Track all combos of residues bridged by solvent.
    typedef std::vector<int> HBlistType;
    HBlistType Donor_;                 ///< Array of hbond donor atoms (D0, H0, D1, H1, ...)
    HBlistType Acceptor_;              ///< Array of hbond acceptor atoms (A0, A1, ...)
    HBlistType SolventDonor_;
    HBlistType SolventAcceptor_;
    AtomMask Mask_;
    AtomMask DonorMask_;
    AtomMask AcceptorMask_;
    AtomMask SolventDonorMask_;
    AtomMask SolventAcceptorMask_;
    bool hasDonorMask_;
    bool hasAcceptorMask_;
    bool hasSolventDonor_;
    bool hasSolventAcceptor_;
    bool calcSolvent_;
    double acut_;
    double dcut2_;

    DataSet* NumHbonds_;
    DataSet* NumSolvent_;
    DataSet* NumBridge_; 
    /// Return true if the first hbond has more frames than the second.
    /// If both have the same # of frames, pick something arbitrary just to
    /// be sure that we have a well-defined ordering (otherwise we could get
    /// spurious test failures) -- Order equivalent frames based on atom number
    /// of Acceptor
    struct hbond_cmp {
      inline bool operator()(HbondType first, HbondType second) const {
        if (first.Frames > second.Frames)
          return true;
        else if (first.Frames < second.Frames)
          return false;
        else
          return (first.A < second.A);
      }
    };

    void SearchAcceptor(HBlistType&,AtomMask&,bool);
    void SearchDonor(HBlistType&,AtomMask&,bool);
    inline int AtomsAreHbonded(int, int, int, int,bool);
};
#endif
