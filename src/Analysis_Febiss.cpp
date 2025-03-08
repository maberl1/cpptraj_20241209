//
// Created by lum on 27.09.23.
//

#include "Analysis_Febiss.h"
#include "CpptrajStdio.h"
#include <iostream>
#include "ProgressBar.h"
#include <memory> //needed by make_shared and unique_ptr

/** This code is based on the GIGist code: https://github.com/liedllab/gigist.git (27 Sept 2023) as well as on the GIST code of cpptraj: https://github.com/Amber-MD/cpptraj.git
 * Action_GIGIST.cpp: https://github.com/liedllab/gigist/blob/9be781be4099f559a2413a553d3464f48ee32a57/Action_GIGIST.cpp (27 Sept 2023)
 * Action_GIGIST.h: https://github.com/liedllab/gigist/blob/9be781be4099f559a2413a553d3464f48ee32a57/Action_GIGIST.h (27 Sept 2023)
 * Action_GIST.cpp: https://github.com/Amber-MD/cpptraj/blob/a60fbeedf3919cc94717cb31299e60eaed800038/src/Action_GIST.cpp (27 Sept 2023)
 * Action_GIST.h: https://github.com/Amber-MD/cpptraj/blob/a60fbeedf3919cc94717cb31299e60eaed800038/src/Action_GIST.h (27 Sept 2023)
*/


/**
 * Standard constructor
 * 
 * Adopted from Action_GIGIST.cpp and edited
 */
Analysis_Febiss::Analysis_Febiss() :
dict_(DataDictionary()),
febissSolventfile_(nullptr)
{}

Analysis_Febiss::Analysis_Febiss(DataSet_3D* population, DataSet_3D* tst, DataSet_3D* tso, DataSet_3D* esw, DataSet_3D* eww, CpptrajFile* out, double refdens, int nframes, unsigned int nsolvent) :
dict_(DataDictionary()),
febissSolventfile_(out)
{
info_.solvent.rho0 = refdens; /** from GIGist::getSystemInfo*/
info_.solvent.numberSolvent = nsolvent;
info_.system.nFrames = nframes;

  result_.resize(25); /** number of dict_ entries*/
  result_.at(dict_.getIndex("population")) = population;
  result_.at(dict_.getIndex("dTStrans_norm")) = tst;
  result_.at(dict_.getIndex("dTSorient_norm")) = tso;
  result_.at(dict_.getIndex("Esw_norm")) = esw;
  result_.at(dict_.getIndex("Eww_norm")) = eww;

# ifdef DEBUG_FEBISS
  mprintf("\t-> Data handling complete!\n");
# endif

  /** ---------------
   * Other variables:
   * ----------------*/
  info_.grid.center = result_.at(dict_.getIndex("population"))->Bin().GridCenter(); /** new: coupled to read in data. initially from GIGist::buildGrid */
  info_.grid.dimensions[0] = result_.at(dict_.getIndex("population"))->NX(); /** new: coupled to read in data. initially from GIGist::buildGrid*/
  info_.grid.dimensions[1] = result_.at(dict_.getIndex("population"))->NY(); /** new: coupled to read in data. initially from GIGist::buildGrid*/
  info_.grid.dimensions[2] = result_.at(dict_.getIndex("population"))->NZ(); /** new: coupled to read in data. initially from GIGist::buildGrid*/
  info_.grid.voxelSize = result_.at(dict_.getIndex("population"))->Bin().DX(); /** new: coupled to read in data. initially from GIGist::buildGrid*/
  info_.grid.voxelVolume = result_.at(dict_.getIndex("population"))->Bin().VoxelVolume(); /** new: coupled to read in data. initially from GIGist::buildGrid*/
  info_.grid.nVoxels = info_.grid.dimensions[0]*info_.grid.dimensions[1]*info_.grid.dimensions[2]; /** new: coupled to read in data. initially from GIGist::buildGrid*/
  info_.grid.start.SetVec(result_.at(dict_.getIndex("population"))->Bin().Corner(0,0,0)[0], result_.at(dict_.getIndex("population"))->Bin().Corner(0,0,0)[1], result_.at(dict_.getIndex("population"))->Bin().Corner(0,0,0)[2]); /** TODO: Where from?*/
  info_.grid.end.SetVec(
      result_.at(dict_.getIndex("population"))->Bin().Corner(info_.grid.dimensions[0],info_.grid.dimensions[1],info_.grid.dimensions[2])[0],
      result_.at(dict_.getIndex("population"))->Bin().Corner(info_.grid.dimensions[0],info_.grid.dimensions[1],info_.grid.dimensions[2])[1],
      result_.at(dict_.getIndex("population"))->Bin().Corner(info_.grid.dimensions[0],info_.grid.dimensions[1],info_.grid.dimensions[2])[2]
      ); /** new: coupled to read in data. initially from GIGist::buildGrid*/

# ifdef DEBUG_FEBISS
  mprintf("\t-> Other variables handling complete!\n");


  mprintf("\n\tSystem informations:\n");
  mprintf("\t\tinfo_.solvent.rho0 = %8.3f\n",info_.solvent.rho0);
  mprintf("\t\tinfo_.solvent.numberSolvent = %i\n",info_.solvent.numberSolvent);
  mprintf("\t\tinfo_.system.nFrames = %i\n",info_.system.nFrames);
  mprintf("\t\tinfo_.grid.center = %8.3f,%8.3f,%8.3f\n",info_.grid.center[0],info_.grid.center[1],info_.grid.center[2]);
  mprintf("\t\tinfo_.grid.dimensions = %i,%i,%i\n",info_.grid.dimensions[0],info_.grid.dimensions[1],info_.grid.dimensions[2]);
  mprintf("\t\tinfo_.grid.voxelSize = %8.3f\n",info_.grid.voxelSize);
  mprintf("\t\tinfo_.grid.voxelVolume = %8.3f\n",info_.grid.voxelVolume);
  mprintf("\t\tinfo_.grid.nVoxels = %i\n",info_.grid.nVoxels);
  mprintf("\t\tinfo_.grid.start = %8.3f,%8.3f,%8.3f\n",info_.grid.start[0],info_.grid.start[1],info_.grid.start[2]);
  mprintf("\t\tinfo_.grid.end = %8.3f,%8.3f,%8.3f\n\n",info_.grid.end[0],info_.grid.end[1],info_.grid.end[2]);
# endif
}


void Analysis_Febiss::Help() const {
mprintf("     Usage:\n"
     " refdens (default: 0.0329 (TIP3P water))"
     " solvnum (total number of solvent molecules)"
     " nframes (total number of frames considered in the previous GIST analysis)"


            "When using FEBISS please cite:"
            "# Steiner, M.; Holzknecht, T.; Schauperl, M.; Podewitz, M. Quantum Chemical Microsolvation by Automated Water Placement. Molecules 2021, 26 (6), 1793. https://doi.org/10.3390/molecules26061793."
            "");
}

// Analysis_Febiss::Setup()
Analysis::RetType Analysis_Febiss::Setup(ArgList& analyzeArgs, ActionInit& setup, int debugIn)
{
/** Arglist handling*/
info_.solvent.rho0 = analyzeArgs.getKeyDouble("refdens", 0.0329); /** from GIGist::getSystemInfo*/
info_.solvent.numberSolvent = analyzeArgs.getKeyInt("solvnum",-1);
info_.system.nFrames = analyzeArgs.getKeyInt("nframes",-1);
std::string febissfile = analyzeArgs.GetStringKey("out", "febiss.dat"); /** new: adopted from Action_GIST.cpp*/

# ifdef DEBUG_FEBISS
mprintf("\t-> Arglist handling complete!\n");
# endif

/** Data handling*/
result_.resize(25); /** number of dict_ entries*/
result_.at(dict_.getIndex("population")) = (DataSet_3D*)setup.DSL().GetDataSet("gist-population.dx"); /** //renamed to "gist-.." LM20231122 //Without (DataSet_3D*) not working, see Analysis_CrankShaft for example //TODO: Catch needed! //renamed from _norm to _dens LM20230813 */
result_.at(dict_.getIndex("dTStrans_norm")) = (DataSet_3D*)setup.DSL().GetDataSet("gist-dTStrans_norm.dx"); /** //renamed to "gist-". Also changed to "-dens" instead of "_dens" LM20231122 //renamed from _norm to _dens LM20230813*/
result_.at(dict_.getIndex("dTSorient_norm")) = (DataSet_3D*)setup.DSL().GetDataSet("gist-dTSorient_norm.dx"); /** //renamed to "gist-". Also changed to "-dens" instead of "_dens" LM20231122  LM20231122 //renamed from _norm to _dens LM20230813*/
result_.at(dict_.getIndex("Esw_norm")) = (DataSet_3D*)setup.DSL().GetDataSet("gist-Esw_norm.dx"); /** //renamed to "gist-". Also changed to "-dens" instead of "_dens" LM20231122  LM20231122 //renamed from _norm to _dens LM20230813*/
result_.at(dict_.getIndex("Eww_norm")) = (DataSet_3D*)setup.DSL().GetDataSet("gist-Eww_norm.dx"); /** //renamed to "gist-". Also changed to "-dens" instead of "_dens" LM20231122  LM20231122 //renamed from _norm to _dens LM20230813*/

# ifdef DEBUG_FEBISS
mprintf("\t-> Data handling complete!\n");
# endif

/** ---------------
 * Other variables:
 * ----------------*/
info_.grid.center = result_.at(dict_.getIndex("population"))->Bin().GridCenter(); /** new: coupled to read in data. initially from GIGist::buildGrid */
info_.grid.dimensions[0] = result_.at(dict_.getIndex("population"))->NX(); /** new: coupled to read in data. initially from GIGist::buildGrid*/
info_.grid.dimensions[1] = result_.at(dict_.getIndex("population"))->NY(); /** new: coupled to read in data. initially from GIGist::buildGrid*/
info_.grid.dimensions[2] = result_.at(dict_.getIndex("population"))->NZ(); /** new: coupled to read in data. initially from GIGist::buildGrid*/
info_.grid.voxelSize = result_.at(dict_.getIndex("population"))->Bin().DX(); /** new: coupled to read in data. initially from GIGist::buildGrid*/
info_.grid.voxelVolume = result_.at(dict_.getIndex("population"))->Bin().VoxelVolume(); /** new: coupled to read in data. initially from GIGist::buildGrid*/
info_.grid.nVoxels = info_.grid.dimensions[0]*info_.grid.dimensions[1]*info_.grid.dimensions[2]; /** new: coupled to read in data. initially from GIGist::buildGrid*/
info_.grid.start.SetVec(result_.at(dict_.getIndex("population"))->Bin().Corner(0,0,0)[0], result_.at(dict_.getIndex("population"))->Bin().Corner(0,0,0)[1], result_.at(dict_.getIndex("population"))->Bin().Corner(0,0,0)[2]); /** TODO: Where from?*/
info_.grid.end.SetVec(
    result_.at(dict_.getIndex("population"))->Bin().Corner(info_.grid.dimensions[0],info_.grid.dimensions[1],info_.grid.dimensions[2])[0],
    result_.at(dict_.getIndex("population"))->Bin().Corner(info_.grid.dimensions[0],info_.grid.dimensions[1],info_.grid.dimensions[2])[1],
    result_.at(dict_.getIndex("population"))->Bin().Corner(info_.grid.dimensions[0],info_.grid.dimensions[1],info_.grid.dimensions[2])[2]
    ); /** new: coupled to read in data. initially from GIGist::buildGrid*/

# ifdef DEBUG_FEBISS
mprintf("\t-> Other variables handling complete!\n");
# endif

/** ------------- 
 * File handling:
 * -------------*/
this->febissSolventfile_ = setup.DFL().AddCpptrajFile(febissfile, "Febiss output"); /** new LM20231213: name for febiss.dat can now be defined. new LM20231122: changed filename to febiss.dat. replaced "createDatasets(argList, actionInit);" in Action_GIGIST.cpp by this->febisssolventfile_ =  ...*/

# ifdef DEBUG_FEBISS
mprintf("\t-> File handling complete!\n");
# endif

//debug
# ifdef DEBUG_FEBISS
  mprintf("\n\tSystem informations:\n");
  mprintf("\t\tinfo_.solvent.rho0 = %8.3f\n",info_.solvent.rho0); 
  mprintf("\t\tinfo_.solvent.numberSolvent = %i\n",info_.solvent.numberSolvent); 
  mprintf("\t\tinfo_.system.nFrames = %i\n",info_.system.nFrames); 
  mprintf("\t\tinfo_.grid.center = %8.3f,%8.3f,%8.3f\n",info_.grid.center[0],info_.grid.center[1],info_.grid.center[2]); 
  mprintf("\t\tinfo_.grid.dimensions = %i,%i,%i\n",info_.grid.dimensions[0],info_.grid.dimensions[1],info_.grid.dimensions[2]); 
  mprintf("\t\tinfo_.grid.voxelSize = %8.3f\n",info_.grid.voxelSize); 
  mprintf("\t\tinfo_.grid.voxelVolume = %8.3f\n",info_.grid.voxelVolume); 
  mprintf("\t\tinfo_.grid.nVoxels = %i\n",info_.grid.nVoxels); 
  mprintf("\t\tinfo_.grid.start = %i,%i,%i\n",info_.grid.start[0],info_.grid.start[1],info_.grid.start[2]); 
  mprintf("\t\tinfo_.grid.end = %i,%i,%i\n",info_.grid.end[0],info_.grid.end[1],info_.grid.end[2]); 
# endif

return Analysis::OK;
}

Analysis::RetType Analysis_Febiss::Analyze() {
placeFebissSolvents();

return Analysis::OK;
}

/**
 * @brief main function for FEBISS placement
 *
 * Adopted from Action_GIGIST.cpp and edited.
 */
void Analysis_Febiss::placeFebissSolvents(void) {
  determineGridShells();

  #ifdef DEBUG_FEBISS
  mprintf("\n\n");
  #endif

  /* calculate delta G and read density data */
  std::vector<double> deltaG;
  std::vector<double> pop;
  for (int voxel = 0; voxel < info_.grid.nVoxels; ++voxel) {
    double dTSt = result_.at(dict_.getIndex("dTStrans_norm"))->operator[](voxel); /** renamed _norm to _dens LM20230813*/
    double dTSo = result_.at(dict_.getIndex("dTSorient_norm"))->operator[](voxel); /** renamed _norm to _dens LM20230813*/
    double esw = result_.at(dict_.getIndex("Esw_norm"))->operator[](voxel); /** renamed _norm to _dens LM20230813*/
    double eww = result_.at(dict_.getIndex("Eww_norm"))->operator[](voxel); /** renamed _norm to _dens LM20230813*/
    double value = esw + eww - dTSo - dTSt;
    deltaG.push_back(value);
    pop.push_back(result_.at(dict_.getIndex("population"))->operator[](voxel)/info_.system.nFrames); /** LM2031206: Added division again, since population is not relpop in GIGist. LM20231130: took back division by rho0, volume, ...//  edited*/
    }
  /* Place solvents to recover 95% of the original density */
  int solventToPosition = static_cast<int>(round(info_.solvent.numberSolvent * 0.95 / 3)); /** edited*/
  mprintf("Placing %d FEBISS solvents\n", solventToPosition);
  ProgressBar progBarFebiss(solventToPosition);
  /* cycle to position all solvents */
  for (int i = 0; i < solventToPosition; ++i) {
    progBarFebiss.Update(i);
    double densityValueOld = 0.0;
    /* get data of current highest density voxel */
    std::vector<double>::iterator maxDensityIterator;
    maxDensityIterator = std::max_element(pop.begin(), pop.end());
    double densityValue = *(maxDensityIterator);
    int max_dens_voxel = maxDensityIterator - pop.begin();
    Vec3 voxelCoords = coordsFromIndex(max_dens_voxel);

    #ifdef DEBUG_FEBISS
    mprintf("\n\ndensityValue %8.3f at %i (%8.3f, %8.3f, %8.3f)\n", densityValue, max_dens_voxel, voxelCoords[0], voxelCoords[1], voxelCoords[2]);
    #endif

    /** Skipping part of Action_GIGIST.cpp that finds the correct rotation of the whole solvent molecule*/

    /* increase included shells until enough density to subtract */
    int shellNum = 0;
    int maxShellNum = shellcontainerKeys_.size() - 1;
    /* not enough density and not reached limit */
    while (densityValue < 1 && /** edited*/
           shellNum < maxShellNum
    ) {
      densityValueOld = densityValue;
      ++shellNum; /** INFO: stepwise increasing of shellnumber, i.e. going from the center of the grid to the border. LM231001*/
      /* new density by having additional solventshell */
      densityValue = addSolventShell(densityValue, pop, max_dens_voxel, shellNum); /** edited*/
    }

    #ifdef DEBUG_FEBISS
    mprintf("\nIn total %i shells added for voxel %i (%8.3f, %8.3f, %8.3f)\nBuilding up weighted deltaG:\n", maxShellNum, max_dens_voxel, voxelCoords[0], voxelCoords[1], voxelCoords[2]);
    #endif

    /* determine density weighted delta G with now reached density */
    double weightedDeltaG = assignDensityWeightedDeltaG(
        max_dens_voxel, shellNum, densityValue, densityValueOld, pop, deltaG);
    /* write new solvent center to pdb */
    writeout(max_dens_voxel, voxelCoords, weightedDeltaG); /** edited*/
    /* subtract density in included shells */
    subtractSolvent(pop, max_dens_voxel, shellNum, densityValue, densityValueOld);
    // cycle of placed solvent molecules

  }
}

/**
 * @brief calculates distance of all voxels to the grid center and groups
 *        them intp shells of identical distances
 *
 * a map stores lists of indices with their identical squared distance as key
 * the indices are stores as difference to the center index to be applicable for
 * all voxels later without knowing the center index
 * a list contains all keys in ascending order to systematically grow included
 * shells later in the algorithm
 * 
 * Adopted from Action_GIGIST.cpp
 */
void Analysis_Febiss::determineGridShells(void) { /** INFO:This is used in placeFebissSolvents */
    /* determine center index */
    size_t centeri, centerj, centerk;
    result_.at(dict_.getIndex("population"))->
            Bin().Calc(info_.grid.center[0], info_.grid.center[1], info_.grid.center[2], centeri, centerj, centerk);
    int centerIndex = result_.at(dict_.getIndex("population"))
            ->CalcIndex(centeri, centerj, centerk);
    /* do not use center_ because it does not align with a voxel but lies between voxels */
    /* however the first shell must be solely the voxel itself -> use coords from center voxel */
    Vec3 centerCoords = coordsFromIndex(centerIndex);
    for (int vox = 0; vox < info_.grid.nVoxels; ++vox) {
        /* determine squared distance */
        Vec3 coords = coordsFromIndex(vox);
        Vec3 difference = coords - centerCoords;
        double distSquared = difference[0] * difference[0] +
                             difference[1] * difference[1] +
                             difference[2] * difference[2];
        /** find function of map returns last memory address of map if not found
        if is entered if distance already present as key in map -> can be added */
        if (shellcontainer_.find(distSquared) != shellcontainer_.end()) {
            shellcontainer_[distSquared].push_back(vox-centerIndex);
        } else {
            /** create new entry in map */
            std::vector<int> indexDifference; /** vector is being reassigned every for-iteration but is needed because the shellcontainer_ is defined as <double,vector<int>>*/
            indexDifference.push_back(vox-centerIndex);
            shellcontainer_.insert(std::make_pair(distSquared, indexDifference));
        }
    }
    /** create list to store ascending keys */
    shellcontainerKeys_.reserve(shellcontainer_.size());
    std::map<double, std::vector<int>>::iterator it = shellcontainer_.begin();
    while(it != shellcontainer_.end()) {
        shellcontainerKeys_.push_back(it->first);
        #ifdef DEBUG_FEBISS
        mprintf("\n\tshellcontainerKey %8.3f: %i values", it->first, it->second.size());
        #endif
        it++;
    }
}

/**
 * @brief own utility function to get coords from grid index
 *
 * @argument index The index in the GIST grid
 * 
 * Adopted from Action_GIGIST.cpp.
 * @return Vec3 coords at the voxel
 */
Vec3 Analysis_Febiss::coordsFromIndex(const int index) {
  size_t i, j, k;
  result_.at(dict_.getIndex("population"))->ReverseIndex(index, i, j, k);
  Vec3 coords = info_.grid.start;
  /* the + 0.5 * size is necessary because of cpptraj's interprets start as corner of grid */
  /* and voxel coordinates are given for the center of the voxel -> hence the shift of half spacing */
  coords[0] += (i + 0.5) * info_.grid.voxelSize;
  coords[1] += (j + 0.5) * info_.grid.voxelSize;
  coords[2] += (k + 0.5) * info_.grid.voxelSize;
  return coords;
}

/**
 * @brief adds density of additional solvent shell to the densityValue
 *
 * @argument densityValue The value after the last shell was added
 * @argument pop The list of the population values
 * @argument index The index of the GIST grid where the central atom will be placed
 * @argument shellNum The number of the new shells to be added
 * 
 * Adopted from Action_GIGIST.cpp::addWaterShell and edited.
 * @return double density value with the new solvent shell
 */
double Analysis_Febiss::addSolventShell(double& densityValue, const std::vector<double>& pop, const int index, const int shellNum) {
    /* get shell from map */
    /** INFO: 
     * shellcontainerKeys_ is a vector that contains the distances of the shells to the center voxel going from the center to the outside. 
     * shellcontainer_ is a map that stores the voxels that belong to shell (i.e. voxels with the same distance to the center voxel). 
     *      The voxel information is given relative to the center voxel as index difference (shellvoxel_index - centervoxel_index). 
     *      This allows to also define shells around voxels other than the center as is done here.*/
    auto newShell = std::make_shared<std::vector<int>>(
            shellcontainer_[shellcontainerKeys_[shellNum]]
    );

    for (unsigned int i = 0; i < (*newShell).size(); ++i) { /** iterates through the voxels that belong to that shell*/
        int tmpIndex = index + (*newShell)[i];
        if (0 < tmpIndex && tmpIndex < static_cast<int>(info_.grid.nVoxels))
            densityValue += pop[tmpIndex];
    }

    #ifdef DEBUG_FEBISS
    mprintf("\tAdded shell %i with distance %8.3f. Density value now: %8.3f\n",shellNum, sqrt(shellcontainerKeys_[shellNum]), densityValue);
    #endif

    return densityValue;
}

/**
 * @brief weights Delta G of GIST with the solvent density that was subtracted
 *        to place the solvent molecule
 *
 * @argument index The index where the central atom is placed
 * @argument shellNum The number of shells around the placed central atom to add to the density
 * @argument densityValue The value after the last shell was added
 * @argument densityValueOld The value before the last shell was added
 * @argument pop The list of the population values
 * @argument deltaG The list of all DeltaG values
 *
 * 
 * Adopted from Action_GIGIST.cpp and edited
 * @return double density weighted Delta G
 */
double Analysis_Febiss::assignDensityWeightedDeltaG(
  int index,
  int shellNum,
  double densityValue,
  double densityValueOld,
  const std::vector<double>& pop,
  const std::vector<double>& deltaG
) {
  double value = 0.0; // value to be returned
  /* cycle through all shells but the last one */
  for (int i = 0; i < shellNum; ++i) {
    /* current shell */
    auto shell = std::make_shared<std::vector<int>>(
        shellcontainer_[shellcontainerKeys_[i]]);
    /* cycle through current shell */
    for (unsigned int j = 0; j < (*shell).size(); ++j) {
      /* get index and check if inside the grid */
      int tmpIndex = index + (*shell)[j];
      if (0 < tmpIndex && tmpIndex < static_cast<int>(info_.grid.nVoxels))
        /* add density weighted delta G to value */
        value += (pop[tmpIndex]/(info_.grid.voxelVolume*info_.solvent.rho0)) * deltaG[tmpIndex];

    }

    #ifdef DEBUG_FEBISS
    mprintf("%i: %8.3f   ", i, value * info_.grid.voxelVolume * info_.solvent.rho0);
    #endif

  }


  double last_shell = densityValue - densityValueOld; // density of last shell
  /* Get percentage of how much of the last shell shall be accounted for */
  double percentage = 1.0;
  if (last_shell != 0.0)
    percentage -= (densityValue - 1) / last_shell; /** edited*/
  /* identical to above but only last shell and percentage */
  auto outerShell = std::make_shared<std::vector<int>>(
      shellcontainer_[shellcontainerKeys_[shellNum]]);
  for (unsigned int i = 0; i < (*outerShell).size(); ++i) {
    int tmpIndex = index + (*outerShell)[i];
    if (0 < tmpIndex && tmpIndex < static_cast<int>(info_.grid.nVoxels))
      value += percentage * pop[tmpIndex]/(info_.grid.voxelVolume*info_.solvent.rho0) * deltaG[tmpIndex];
  }

    #ifdef DEBUG_FEBISS
    mprintf("%i (outer shell): %8.3f (percentage = %8.3f)   ", shellNum, value * info_.grid.voxelVolume * info_.solvent.rho0, percentage);
    #endif

  return value * info_.grid.voxelVolume * info_.solvent.rho0; /** edited*/
}

/**
 * @brief writes out the voxel, its coords and assigned deltaG value
 *
 * @argument atomNumber The running index in the pdb file
 * @argument voxelCoords The coordinates for the central atom
 * @argument deltaG density weighted Delta G to be included as b-factor
 * 
 * Adopted from Action_GIGIST.cpp (Action_Febiss::writeFebissPdb) and edited.
 */
void Analysis_Febiss::writeout
        (
                const int max_dens_voxel,
                const Vec3& voxelCoords,
                const double deltaG
        ) {

    febissSolventfile_->Printf( /** change output to: tuple(voxelnumber, voxelCoords[0], voxelCoords[1], voxelCoords[2], weighted deltaG)*/
            "%d %8.3f%8.3f%8.3f%7.2f\n", 
            max_dens_voxel,
            voxelCoords[0],
            voxelCoords[1],
            voxelCoords[2],
            deltaG
    );
}

/**
 * @brief subtract density from all voxels that belonged to the included shells
 *
 * @argument pop The list of the population values
 * @argument index The index where the central atom is placed
 * @argument shellNum The number of shells around the placed central atom that were included
 * @argument densityValue The value after the last shell was added
 * @argument densityValueOld The value before the last shell was added
 * 
 * Adopted from Action_Febiss.cpp and edited.
 */
void Analysis_Febiss::subtractSolvent
        (
                std::vector<double>& pop,
                int index,
                int shellNum,
                double densityValue,
                double densityValueOld
        )
{
    /* cycle through all but the last shell */
    for (int i = 0; i < shellNum; ++i) {
        auto shell = std::make_shared<std::vector<int>>(
                shellcontainer_[shellcontainerKeys_[i]]);
        for (unsigned int j = 0; j < (*shell).size(); ++j) {
            int tmpIndex = index + (*shell)[j];
            if (0 < tmpIndex && tmpIndex < static_cast<int>(info_.grid.nVoxels))
                /* remove all population from the voxel in the GIST grid */
                pop[tmpIndex] = 0.0;
        }
    }
    /* since density of one solvent is overshot, the density must not be deleted
     * completely in the last shell, first determine percentage */
    double last_shell = densityValue - densityValueOld; // density in last shell
    double percentage = 1.0;
    if (last_shell != 0.0)
        percentage -= (densityValue - 1)/ last_shell;
    /* identical to before but only last shell and percentage */
    auto outerShell = std::make_shared<std::vector<int>>(
            shellcontainer_[shellcontainerKeys_[shellNum]]);
    for (unsigned int i = 0; i < (*outerShell).size(); ++i) {
        int tmpIndex = index + (*outerShell)[i];
        if (0 < tmpIndex && tmpIndex < static_cast<int>(info_.grid.nVoxels))
            pop[tmpIndex] -= percentage * pop[tmpIndex];
    }
}