//
// Created by lum on 27.09.23.
//

#ifndef INC_ANALYSIS_FEBISS_H
#define INC_ANALYSIS_FEBISS_H
#include "Analysis.h"
#include "DataSet_3D.h"
#include <map>
#include <array>
#include "ProgressBar.h"
//#include "DataSet_GridFlt.h" /** inherits <algorithm> for std::max_element*/
#include <algorithm>


/** DataDictionary class adopted from Action_GIGist.h: https://github.com/liedllab/gigist/blob/9be781be4099f559a2413a553d3464f48ee32a57/Action_GIGIST.h (27.Sept.2023) */

/**
 * Data Dictionary helper class.
 * The different atoms of the solvent can be added to the
 * dictionary via the add command.
 */
class DataDictionary {
private:
  std::vector<std::string> names;
public:
  /**
   * Constructor creates the initial data points, which will always be the same
   * for each GIST run.
   */
  DataDictionary() {
    this->names.push_back("population");
    this->names.push_back("dTStrans_norm");
    this->names.push_back("dTStrans_dens");
    this->names.push_back("dTSorient_norm");
    this->names.push_back("dTSorient_dens");
    this->names.push_back("dTSsix_norm");
    this->names.push_back("dTSsix_dens");
    this->names.push_back("Eww");
    this->names.push_back("Eww_norm");
    this->names.push_back("Eww_dens");
    this->names.push_back("Esw");
    this->names.push_back("Esw_norm");
    this->names.push_back("Esw_dens");
    this->names.push_back("dipole_x");
    this->names.push_back("dipole_y");
    this->names.push_back("dipole_z");
    this->names.push_back("dipole_xtemp");
    this->names.push_back("dipole_ytemp");
    this->names.push_back("dipole_ztemp");
    this->names.push_back("dipole_g");
    this->names.push_back("order");
    this->names.push_back("order_norm");
    this->names.push_back("neighbour");
    this->names.push_back("neighbour_dens");
    this->names.push_back("neighbour_norm");
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
        struct Solute {
            int numberSoluteAtoms = 0; /** built up in DoAction -> writeOutSolute*/
        } solute;
        struct Solvent {
            double rho0 = 0.0; /** defined in Setup.*/
            std::string centralatom = "O"; /** new: label of central atom to be placed. default: O for use with water*/
            int solventAtomCounter_; /** adopted from Action_GIGIST.h. There it is a standalone std::vector<int> variable.*/
            //int rigidAtomsN = 3; /** new: how many rigidatoms shall be defined. TODO: make it constant 3 if needed at all*/
            int numberSolvent = 0; /** defined in Setup.*/
            //int numberSolventAtoms = 0; /** defined in Init.*/
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

  std::vector<DataSet_3D*> result_; /** Adopted from Action_GIGIST.h. This variable is where GIST data is loaded into.*/
  CpptrajFile *febissSolventfile_; /** Adopted from Action_GIGIST.h. This file contains information on the energy and the position of the solvent molecules.*/
  std::map<double, std::vector<int>> shellcontainer_; /** Adopted from Action_GIGIST.h. */
  std::vector<double> shellcontainerKeys_; /** Adopted from Action_GIGIST.h. */
  DataDictionary dict_; /** Adopted from Action_GIGIST.h. */


  /** -----------
      * Functions
      *----------*/

  Vec3 coordsFromIndex(const int); /** Adopted from Action_GIGIST.h */

  void placeFebissSolvents(void); /** Adopted from Action_GIGIST.h */
  
  void determineGridShells(void); /** Adopted from Action_GIGIST.h */

  double assignDensityWeightedDeltaG(const int, const int, const double, const double, const std::vector<double>&, const std::vector<double>&);

  double addSolventShell(double&, const std::vector<double>&, const int, const int);

  void writeout(const int, const Vec3&, const double); /** Adopted from Action_GIGIST.h and edited.*/

  void subtractSolvent(std::vector<double> &, const int, const int, const double, const double); /** Adopted from Action_GIGIST.h and edited.*/

};
#endif

