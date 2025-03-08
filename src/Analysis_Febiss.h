#ifndef INC_ANALYSIS_FEBISS_H
#define INC_ANALYSIS_FEBISS_H
#include "Analysis.h"
#include "DataSet_3D.h"
#include <map>
#include <array>
#include "ProgressBar.h"
#include <algorithm>


/** DataDictionary class adopted from Action_GIGist.h: https://github.com/liedllab/gigist/blob/9be781be4099f559a2413a553d3464f48ee32a57/Action_GIGIST.h (27.Sept.2023) */

/**
 * Data Dictionary helper class.
 */
class DataDictionary {
private:
  std::vector<std::string> names;

public:

  DataDictionary() {
    this->names.push_back("population");
    this->names.push_back("dTStrans_norm");
    this->names.push_back("dTSorient_norm");
    this->names.push_back("dTSsix_norm");
    this->names.push_back("Eww_norm");
    this->names.push_back("Esw_norm");
  }

  /**
   * Calculate the number of data sets.
   * @return: The size of the dictionary.
   */
  unsigned int size( void ) {
    return this->names.size();
  }

  /**
   * Returns the index of a given name, -1, if it is not present.
   * @argument testString: The string to be checked.
   * @return: The index of the testString, or -1 if testString is
   *           not in the data sets.
   */
  int getIndex(std::string testString) {
    for (unsigned int i = 0; i < this->names.size(); ++i) {
      if (testString.compare(this->names.at(i)) == 0) {
        return i;
      }
    }
    return -1;
  }
#if FALSE
  /**
   * Checks whether the dictionary already contains a specific data set.
   * @argument testString: The string to be checked.
   * @return: True if a data set with this name is already there, false otherwise.
   */
  bool contains(std::string testString) {
    return (this->getIndex(testString) != -1);
  }

  /**
   * Get the name of the data set at a given index.
   * @argument idx: The index at which the data set is present.
   * @return: The name of the data set.
   */
  std::string getElement(int idx) {
    return this->names.at(idx);
  }

  /**
   * Add a data set name to the dictionary.
   * @argument name: The name of the new data set.
   */
  void add(std::string name) {
    this->names.push_back(name);
  }
#endif
};


class Analysis_Febiss : public Analysis {
  public:
    Analysis_Febiss();
    Analysis_Febiss(DataSet_3D*, DataSet_3D*, DataSet_3D*, DataSet_3D*, DataSet_3D*, CpptrajFile*, double, int, unsigned int);
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Febiss(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:

  /** ------------
      * Variables
      * ----------*/

  /** 'struct Info' adopted from GIGIST and edited*/
  // Necessary information for the computation
    struct Info {
        struct System {
            int nFrames = 0; /** defined in Setup.*/
        } system;
        struct Solvent {
            double rho0 = 0.0; /** defined in Setup.*/
            int numberSolvent = 0; /** defined in Setup.*/
        } solvent;
        struct Grid {
            double voxelSize = 0.0; /** defined in Setup.*/
            double voxelVolume = 0.0; /** defined in Setup.*/
            std::array<int, 3> dimensions; /** defined in Setup.*/
            int nVoxels = 0; /** defined in Setup.*/
            Vec3 center; /** defined in Setup.*/
            Vec3 start; /** defined in Setup.*/
            Vec3 end; /** defined in Setup.*/
        } grid;
    } info_;

  std::vector<DataSet_3D*> result_; /** This variable is where GIST data is loaded into.*/
  CpptrajFile *febissSolventfile_; /** This file contains information on the energy and the position of the solvent molecules.*/
  std::map<double, std::vector<int>> shellcontainer_; 
  std::vector<double> shellcontainerKeys_; 
  DataDictionary dict_; 


  /** -----------
      * Functions
      *----------*/

  Vec3 coordsFromIndex(const int); 

  void placeFebissSolvents(void); 
  
  void determineGridShells(void); 

  double assignDensityWeightedDeltaG(const int, const int, const double, const double, const std::vector<double>&, const std::vector<double>&);

  double addSolventShell(double&, const std::vector<double>&, const int, const int);

  void writeout(const int, const Vec3&, const double); 

  void subtractSolvent(std::vector<double> &, const int, const int, const double, const double); 

};
#endif

