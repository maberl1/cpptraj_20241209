#include <cmath>
#include <cstdio>
#include "Constants.h"
#include "Matrix_3x3.h"

// MomentOfInertia()
static void MomentOfInertia(int natom, double *X_, double* Mass_, double* pmom)
{
  double IVEC[9];
  double eigvec[9];
  double e2[3];
  // Center of mass
  double cx = 0.0;
  double cy = 0.0;
  double cz = 0.0;
  double sumMass = 0.0;
  double* mass = Mass_;
  int natom3 = natom * 3;
  for (int i = 0; i < natom3; i+=3) {
    sumMass += (*mass);
    cx += ( X_[i  ] * (*mass) );
    cy += ( X_[i+1] * (*mass) );
    cz += ( X_[i+2] * (*mass) );
    ++mass;
  }
  cx /= sumMass; 
  cy /= sumMass; 
  cz /= sumMass;

  // Moment of inertia
  double xx = 0.0; 
  double yy = 0.0; 
  double zz = 0.0; 
  double xy = 0.0; 
  double xz = 0.0; 
  double yz = 0.0; 
  mass = Mass_;
  for (int i = 0; i < natom3; i+=3) {
    double dx = X_[i  ] - cx;
    double dy = X_[i+1] - cy;
    double dz = X_[i+2] - cz;

    xx += *mass * ( dy * dy + dz * dz );
    yy += *mass * ( dx * dx + dz * dz );
    zz += *mass * ( dx * dx + dy * dy );
    xy -= *mass * dx * dy;
    xz -= *mass * dx * dz;
    yz -= *(mass++) * dy * dz;
  }
  IVEC[0] = xx;
  IVEC[1] = xy;
  IVEC[2] = xz;
  IVEC[3] = xy;
  IVEC[4] = yy;
  IVEC[5] = yz;
  IVEC[6] = xz;
  IVEC[7] = yz;
  IVEC[8] = zz;

  Matrix_3x3 TEMP( IVEC );
  // NOTE: Diagonalize sorts evals/evecs in descending order, but
  //       thermo() expects ascending.
  TEMP.Diagonalize_Sort( eigvec, e2 );
  pmom[0] = e2[2];
  pmom[1] = e2[1];
  pmom[2] = e2[0];
  /*int info = 0;
  char jobz = 'V';
  char uplo = 'U';
  int nvec = 3;
  dspev_(jobz, uplo, nvec, IVEC, pmom, eigvec, nvec, e2, info);*/
}

// thermo() 
/** Given the structure of a molecule and its normal mode vibrational
  * frequencies this routine uses standard statistical mechanical
  * formulas for an ideal gas (in the canonical ensemble, see,
  * for example, d. a. mcquarrie, "statistical thermodynamics",
  * harper & row, new york, 1973, chapters 5, 6, and 8) to compute
  * the entropy, heat capacity, and internal energy.

  * The si system of units is used internally. Conversion to units
  * more familiar to most chemists is made for output.
  *
  * \param crd     coordinates in Angstroms
  * \param amass   atomic weights, in amu.
  * \param freq    vibrational frequencies, in cm**-1 and in ascending order
  * \param vtemp   vibrational temperatures, in kelvin.
  * \param evibn   contribution to e from the vibration n.
  * \param cvibn   contribution to cv from the vibration n.
  * \param svibn   contribution to s from the vibration n.
  * \param temp    temperature
  * \param patm    pressure, in atmospheres
*/
void thermo( int natoms, int nvecs, int ilevel, double* crd, double* amass,
             double* freq, double* vtemp, double* evibn, double* cvibn,
             double* svibn, double temp, double patm)
{
  // pmom    principal moments of inertia, in amu-bohr**2 and in ascending order.
  double pmom[10], rtemp, rtemp1, rtemp2, rtemp3;

  // ----- Constants -------------------
  const double thresh = 900.0;        // vibrational frequency threshold
  const double tokg   = 1.660531e-27; // kilograms per amu.
  const double boltz  = 1.380622e-23; // boltzman constant, in joules per kelvin.
  const double planck = 6.626196e-34; // planck constant, in joule-seconds.
  const double avog   = 6.022169e+23; // avogadro constant, in mol**(-1).
  const double jpcal  = 4.18674e+00;  // joules per calorie.
  const double tomet  = 1.0e-10;      // metres per Angstrom.
  const double hartre = 4.35981e-18;  // joules per hartree.
  const double pstd   = 1.01325e+05;  // Standard pressure in pascals
  // -----------------------------------

  //     compute the gas constant, pi, pi**2, and e.
  //     compute the conversion factors cal per joule and kcal per joule.
  const double gas  = avog * boltz;
  // pi   = four * datan(one)
  const double pipi = PI * PI;
  const double e    = exp(1.0);
  const double tocal  = 1.0 / jpcal;
  const double tokcal = tocal / 1000.0;

  //     print the temperature and pressure.
  printf("\n                    *******************\n");
  printf(  "                    - Thermochemistry -\n");
  printf(  "                    *******************\n\n");
  printf("\n temperature %9.3f kelvin\n pressure    %9.5f atm\n",temp,patm);
  double pressure = pstd * patm;
  double rt = gas * temp;

  //     compute and print the molecular mass in amu, then convert to
  //     kilograms.
  double weight = 0.0;
  for (int iat = 0; iat < natoms; ++iat)
    weight += amass[iat];
  printf(" molecular mass (principal isotopes) %11.5f amu\n", weight);
  weight *= tokg;

  //     trap non-unit multiplicities.
  //     if (multip != 1) {
  //       printf("\n Warning-- assumptions made about the electronic partition function\n");
  //       printf(  "           are not valid for multiplets!\n\n");

  //     compute contributions due to translation:
  //        etran-- internal energy
  //        ctran-- constant v heat capacity
  //        stran-- entropy
  double dum1 = boltz * temp;
  double dum2 = pow(TWOPI, 1.5);
  double arg = pow(dum1, 1.5) / planck;
  arg = (arg / pressure) * (dum1 / planck);
  arg = arg * dum2 * (weight / planck);
  arg = arg * sqrt(weight) * exp(2.5);
  double stran = gas * log(arg);
  double etran = 1.5 * rt;
  double ctran = 1.5 * gas;

  //     Compute contributions due to electronic motion:
  //        It is assumed that the first electronic excitation energy
  //        is much greater than kt and that the ground state has a
  //        degeneracy of one.  Under these conditions the electronic
  //        partition function can be considered to be unity.  The
  //        ground electronic state is taken to be the zero of
  //        electronic energy.

  //     for monatomics print and return.
  if (natoms <= 1){ 
    printf("\n internal energy:   %10.3f joule/mol         %10.3f kcal/mol\n",
           etran, etran * tokcal);
    printf(  " entropy:           %10.3f joule/k-mol       %10.3f cal/k-mol\n",
           stran, stran * tocal);
    printf(  " heat capacity cv:  %10.3f joule/k-mol       %10.3f  cal/k-mol\n",
           ctran, ctran * tocal);
    return;
  }

  //     compute contributions due to rotation.

  //     Compute the principal moments of inertia, get the rotational
  //     symmetry number, see if the molecule is linear, and compute
  //     the rotational temperatures.  Note the imbedded conversion
  //     of the moments to SI units.
  MomentOfInertia( natoms, crd, amass, pmom );
  printf("\n principal moments of inertia (nuclei only) in amu-A**2:\n");
  printf(  "      %12.2f%12.2f%12.2f\n", pmom[0], pmom[1], pmom[2]);
  
  bool linear = false;
  // Symmetry number: only for linear molecules. for others symmetry number is unity
  double sn = 1.0;
  if (natoms <= 2) {
    linear = true;
    if (amass[0] == amass[1]) sn = 2.0;
  }
  printf("\n rotational symmetry number %3.0f\n", sn);

  double con = planck / (boltz*8.0*pipi);
  con = (con / tokg)  *  (planck / (tomet*tomet));
  if (linear) {
    rtemp = con / pmom[2];
    if (rtemp < 0.2) {
      printf("\n Warning-- assumption of classical behavior for rotation\n");
      printf(  " may cause significant error\n");
    }
    printf("\n rotational temperature (kelvin) %12.5f\n", rtemp);                 
  } else {
    rtemp1 = con / pmom[0];
    rtemp2 = con / pmom[1];
    rtemp3 = con / pmom[2];
    if (rtemp1 < 0.2) {
      printf("\n Warning-- assumption of classical behavior for rotation\n");
      printf(  " may cause significant error\n");
    }
    printf("\n rotational temperatures (kelvin) %12.5f%12.5f%12.5f\n", 
           rtemp1, rtemp2, rtemp3);
  }

  //         erot-- rotational contribution to internal energy.
  //         crot-- rotational contribution to cv.
  //         srot-- rotational contribution to entropy.
  double erot, crot, srot;

  if (linear) { 
     erot = rt;
     crot = gas;
     arg  = (temp/rtemp) * (e/sn);
     srot = gas * log(arg);
  } else {
     erot = 1.5 * rt;
     crot = 1.5 * gas;
     arg  = sqrt(PI*e*e*e) / sn;
     double dum  = (temp/rtemp1) * (temp/rtemp2) * (temp/rtemp3);
     arg  = arg * sqrt(dum);
     srot = gas * log(arg);
  }

  //     compute contributions due to vibration.

  //     compute vibrational temperatures and zero point vibrational
  //     energy.  only real frequencies are included in the analysis.

  //     ndof = 3*natoms - 6 - nimag
  //     if (nimag .ne. 0) write(iout,1210) nimag
  //     if (linear) ndof = ndof + 1
  int ndof = nvecs;

  //       (---iff is the first frequency to include in thermo:)
  int iff;
  if (ilevel != 0)
     iff = 0;
  else if (linear)
     iff = 5;
  else
     iff = 6;
  con = planck / boltz;
  double ezpe = 0.0;
  for (int i = 0; i < ndof; ++i) {
     vtemp[i] = freq[i+iff] * con * 3.0e10;
     ezpe    += freq[i+iff] * 3.0e10;
  }
  ezpe = 0.5 * planck * ezpe;
  printf("\n zero point vibrational energy %12.1f (joules/mol) \n",ezpe * avog);
  printf(  "                               %12.5f (kcal/mol)\n",ezpe * tokcal * avog);
  printf(  "                               %12.7f (hartree/particle)\n", ezpe / hartre); 
  //     compute the number of vibrations for which more than 5% of an
  //     assembly of molecules would exist in vibrational excited states.
  //     special printing for these modes is done to allow the user to
  //     easily take internal rotations into account.  the criterion
  //     corresponds roughly to a low frequency of 1.9(10**13) hz, or
  //     625 cm**(-1), or a vibrational temperature of 900 k.

  int lofreq = 0;
  for (int i = 0; i < ndof; ++i)
    if (vtemp[i] < thresh)
      ++lofreq;
  if (lofreq != 0) {
    printf("\n Warning-- %3i vibrations have low frequencies and may represent hindered \n",
           lofreq);
    printf(  "         internal rotations.  The contributions printed below assume that these \n");
    printf(  "         really are vibrations.\n");
  }

  //     compute:
  //        evib-- the vibrational component of the internal energy.
  //        cvib-- the vibrational component of the heat capacity.
  //        svib-- the vibrational component of the entropy.
  double evib = 0.0;
  double cvib = 0.0;
  double svib = 0.0;
  double scont;
  for (int i = 0; i < ndof; ++i) {
     //       compute some common factors.

     double tovt  = vtemp[i] / temp;
     double etovt = exp(tovt);
     double em1   = etovt - 1.0;

     //       compute contributions due to the i'th vibration.

     double econt = tovt  *  (0.5 + 1.0/em1);
     double ccont = etovt *  pow(tovt/em1,2.0);
     double argd = 1.0 - 1.0/etovt;
     if (argd > 1.0e-7) 
        scont = tovt/em1 - log(argd);
     else {
        scont = 0.0;
        printf(" warning: setting vibrational entropy to zero for mode %i with vtemp = %f\n",
               i+1, vtemp[i]);
     }
     //       if (lofreq .ge. i) then
     evibn[i] = econt * rt;
     cvibn[i] = ccont * gas;
     svibn[i] = scont * gas;
     //       end if
     evib += econt;
     cvib += ccont;
     svib += scont;
  } 
  evib *= rt;
  cvib *= gas;
  svib *= gas;

  //     the units are now:
  //         e-- joules/mol
  //         c-- joules/mol-kelvin
  //         s-- joules/mol-kelvin

  double etot = etran + erot + evib;
  double ctot = ctran + crot + cvib;
  double stot = stran + srot + svib;

  //     print the sum of the hartree-fock energy and the thermal energy.

  //     call tread(501,gen,47,1,47,1,0)
  //     esum = gen(32) + etot/avog/hartre
  //     write(iout,1230) esum

  //     convert to the following and print
  //         e-- kcal/mol
  //         c-- cal/mol-kelvin
  //         s-- cal/mol-kelvin
  etran = etran * tokcal;
  ctran = ctran * tocal;
  stran = stran * tocal;
  erot   = erot * tokcal;
  crot   = crot * tocal;
  srot   = srot * tocal;
  evib   = evib * tokcal;
  cvib   = cvib * tocal;
  svib   = svib * tocal;
  etot   = etran + erot + evib;
  ctot   = ctran + crot + cvib;
  stot   = stran + srot + svib;
  for (int i = 0; i < ndof; ++i) {
     evibn[i] *= tokcal;
     cvibn[i] *= tocal;
     svibn[i] *= tocal;
  }

  printf("\n\n           freq.         E                  Cv                 S\n");
  printf(    "          cm**-1      kcal/mol        cal/mol-kelvin    cal/mol-kelvin\n");
  printf(    "--------------------------------------------------------------------------------\n");
  printf(    " Total              %11.3f        %11.3f        %11.3f    \n",etot,ctot,stot);
  printf(    " translational      %11.3f        %11.3f        %11.3f    \n",etran,ctran,stran);
  printf(    " rotational         %11.3f        %11.3f        %11.3f    \n",erot,crot,srot);
  printf(    " vibrational        %11.3f        %11.3f        %11.3f    \n",evib,cvib,svib);

  for (int i = 0; i < iff; ++i) 
    printf(" %5i%10.3f\n", i+1, freq[i]);

  for (int i = 0; i < ndof; ++i) {
    printf(" %5i%10.3f    %11.3f        %11.3f        %11.3f    \n",i+iff+1,
           freq[i+iff], evibn[i], cvibn[i], svibn[i]);
  }

}