# This should be sourced at the top of CPPTRAJ test run scripts.
#
# The purpose of this script is to provide all common functionality required
# by all tests, including environment setup and collecting test results.
# Requires awk, grep, sed, diff (if not AmberTools), rm, make (multiple tests)
#
# ----- Environment variables ------------------------------
# Binary locations
#   CPPTRAJ              : Set to cpptraj binary being tested.
#   AMBPDB               : Set to ambpdb binary being tested.
#   VALGRIND             : Set to 'valgrind' command if memory check requested.
#   CPPTRAJ_DIFF         : Command used to check for test differences.
#   CPPTRAJ_DACDIF       : Set if using 'dacdif' in AmberTools to checkout for test differences.
#   CPPTRAJ_NDIFF        : Set to Nelson H. F. Beebe's ndiff.awk for numerical diff calc.
#   CPPTRAJ_NCDUMP       : Set to ncdump command; needed for NcTest()
#   CPPTRAJ_RM           : Command used to remove files
#   CPPTRAJ_TIME         : Set to the 'time' command if timing requested.
#   CPPTRAJ_NPROC        : nproc binary for counting threads in parallel tests.
# Test output locations
#   CPPTRAJ_TEST_RESULTS : File to record individual test results to.
#   CPPTRAJ_TEST_ERROR   : File to record individual test errors/diffs to.
#   CPPTRAJ_OUTPUT       : File to direct cpptraj STDOUT to.
#   CPPTRAJ_ERROR        : File to direct cpptraj STDERR to.
# Other variables
#   CPPTRAJ_TEST_ROOT    : Test root directory. FIXME make local?
#   CPPTRAJ_TEST_SETUP   : 'yes' if setup is complete.
#   CPPTRAJ_TEST_CLEAN   : If 1, only cleaning tests; do not run them.
#   CPPTRAJ_TEST_OS      : Operating system on which tests are being run. If blank assume linux.
#   N_THREADS            : Number of MPI threads if parallel.
#   OMP_NUM_THREADS      : Max number of OpenMP threads.
#   DO_PARALLEL          : MPI run command (e.g. 'mpirun -n 11')
#   CPPTRAJ_DEBUG        : Can be set to pass global debug flag to cpptraj.
#   DIFFOPTS             : Additional options to pass to CPPTRAJ_DIFF
#   CPPTRAJ_PROFILE      : If 1, end of test profiling with gprof performed.
# Cpptraj binary characteristics
#   CPPTRAJ_ZLIB         : If set CPPTRAJ has zlib support.
#   CPPTRAJ_BZLIB        : If set CPPTRAJ has bzip support.
#   CPPTRAJ_NETCDFLIB    : If set CPPTRAJ has NetCDF support.
#   CPPTRAJ_MPILIB       : If set CPPTRAJ has MPI support.
#   CPPTRAJ_NOMATHLIB    : If set CPPTRAJ was compiled without math libraries.
#   CPPTRAJ_OPENMP       : If set CPPTRAJ has OpenMP support.
#   CPPTRAJ_PNETCDFLIB   : If set CPPTRAJ has parallel NetCDF support.
#   CPPTRAJ_SANDERLIB    : If set CPPTRAJ was compiled with the sander API from AT.
#   CPPTRAJ_CUDA         : If set CPPTRAJ has CUDA support.
#   CPPTRAJ_NO_XDRFILE   : If set CPPTRAJ was compiled without XDR file support.
# ----- Variables that can be set by individual scripts ----
#   TOP                  : Topology file for cpptraj
#   INPUT                : Input file for cpptraj
#   CPPTRAJ_TEST_MODE    : Set to 'master' if executed from CpptrajTest.sh.
# ----- Local setup variables ------------------------------
VGMODE=0                 # Valgrind mode: 0 none, 1 memcheck, 2 helgrind
USE_DACDIF=1             # If 0 do not use dacdif even if in AmberTools
GET_TIMING=0             # If 1 time cpptraj with the CPPTRAJ_TIME binary
SUMMARY=0                # If 1 print summary of CPPTRAJ_TEST_RESULTS only
SHOWERRORS=0             # If 1, print test errors to STDOUT after summary.
SFX=""                   # CPPTRAJ binary suffix
TARGET=""                # Make target if multiple tests being run
STANDALONE=1             # If 0, part of AmberTools. If 1, stand-alone (e.g. from GitHub).
TEST_DIRS=''             # Specific test directories to run.
# ----- Variables local to single test ---------------------
TEST_WORKDIR=''          # Test working directory
NUMCOMPARISONS=0         # Total number of times DoTest has been called this test.
ERRCOUNT=0               # Total number of errors detected by DoTest this test.
WARNCOUNT=0              # Total number of warnings detected by DoTest this test.
PROGCOUNT=0              # Total number of times RunCpptraj has been called this test.
PROGERROR=0              # Total number of program errors this test

# ==============================================================================
# TestHeader() <outfile>
#   Write test header (working directory) to specified file.
TestHeader() {
  echo "**************************************************************" > $1
  echo "TEST: $TEST_WORKDIR" >> $1
}

# OUT() <message>
#   Write message to CPPTRAJ_TEST_RESULTS
OUT() {
  echo "$1" >> $CPPTRAJ_TEST_RESULTS
}

# ERR() <message>
#   Write message to CPPTRAJ_TEST_ERROR. Write header if this is the first error.
ERR() {
  if [ $ERRCOUNT -eq 0 -a -z "$CPPTRAJ_DACDIF" ] ; then
    TestHeader "$CPPTRAJ_TEST_ERROR"
  fi
  echo "$1" >> $CPPTRAJ_TEST_ERROR
}

# OutBoth() <message>
#   Send <message> to CPPTRAJ_TEST_RESULTS and CPPTRAJ_TEST_ERROR
OutBoth() {
  OUT "$1"
  ERR "$1"
}

# FIXME This is a stub and should be removed
CheckTest() {
  CHECKTEST=1
}

# ------------------------------------------------------------------------------
# DoTest() <File1> <File2> [-r <relative err>] [-a <absolute err>] [<arg1>] ... [<argN>]
#   Compare File1 (the 'save' file) to File2 (test output), print an error if
#   they differ. If '-r' or '-a' are specified the test will pass as long as the
#   maximum relative or absolulte errors are below the given value (using Nelson
#   H. F. Beebe's ndiff.awk script. The remaining args can be used to pass
#   options to CPPTRAJ_DIFF.
DoTest() {
  #echo "DEBUG: DoTest $1 $2"
  if [ ! -z "$CPPTRAJ_DACDIF" ] ; then
    # Use dacdif. Use any '-r <X>' or '-a <X>' args found. Ignore the rest.
    DIFFARGS="$1 $2"
    shift # Save file
    shift # Test file
    # Process remaining args
    while [ ! -z "$1" ] ; do
      case "$1" in
        "-r" ) shift ; DIFFARGS=" -r $1 "$DIFFARGS ;;
        "-a" ) shift ; DIFFARGS=" -a $1 "$DIFFARGS ;;
      esac
      shift
    done
    $CPPTRAJ_DACDIF $DIFFARGS
  else
    # Use diff, or ndiff where '-r' or '-a' specified.
    ((NUMCOMPARISONS++))
    DIFFARGS='--strip-trailing-cr'
    NDIFFARGS=""
    # First two arguments are files to compare.
    F1=$1 ; shift
    F2=$1 ; shift
    # Process remaining arguments.
    USE_NDIFF=0
    while [ ! -z "$1" ] ; do
      case "$1" in
        "-r"           ) USE_NDIFF=1; shift; NDIFFARGS="$NDIFFARGS -v RELERR=$1" ;;
        "-a"           ) USE_NDIFF=1; shift; NDIFFARGS="$NDIFFARGS -v ABSERR=$1" ;;
        *              ) DIFFARGS=$DIFFARGS" $1" ;;
      esac
      shift
    done
    if [ ! -f "$F1" ] ; then
      OutBoth "  Save file '$F1' not found."
      ((ERRCOUNT++))
    elif [ ! -f "$F2" ] ; then
      OutBoth "  Test output '$F2' not found."
      ((ERRCOUNT++))
    else
      if [ $USE_NDIFF -eq 0 ] ; then
        $CPPTRAJ_DIFF $DIFFARGS $DIFFOPTS $F1 $F2 > temp.diff 2>&1
      else
        awk -f $CPPTRAJ_NDIFF $NDIFFARGS $F1 $F2 > temp.diff 2>&1
      fi
      if [ -s 'temp.diff' ] ; then
        OutBoth "  $F1 $F2 are different."
        cat temp.diff >> $CPPTRAJ_TEST_ERROR
        ((ERRCOUNT++))
      else
        OUT  "  $F2 OK."
      fi
      $CPPTRAJ_RM temp.diff
    fi
  fi
}

# ------------------------------------------------------------------------------
# NcTest() <1> <2>
#   Compare NetCDF files <1> and <2>. Use CPPTRAJ_NCDUMP to convert to ASCII
#   first, removing ==> line and :programVersion attribute.
NcTest() {
  #echo "DEBUG: NcTest $1 $2"
  if [ -z "$1" -o -z "$2" ] ; then
    echo "Error: NcTest(): One or both files not specified." > /dev/stderr
    exit 1
  fi
  # Save remaining args for DoTest
  F1=$1
  F2=$2
  shift
  shift
  DIFFARGS="nc0.save nc0"
  CALC_NUM_ERR=0
  while [ ! -z "$1" ] ; do
    if [ "$1" = '-r' -o "$1" = '-a' ] ; then
      CALC_NUM_ERR=1
    fi
    DIFFARGS=$DIFFARGS" $1"
    shift
  done
  # Prepare files.
  if [ ! -e "$F1" ] ; then
    if [ -z "$CPPTRAJ_DACDIF" ] ; then
      OutBoth "  Save file '$F1' not found."
    fi
    ((NUMCOMPARISONS++))
    ((ERRCOUNT++))
  elif [ ! -e "$F2" ] ; then
    if [ -z "$CPPTRAJ_DACDIF" ] ; then
      OutBoth "  Test output '$F2' not found."
    fi
    ((NUMCOMPARISONS++))
    ((ERRCOUNT++))
  else
    if [ $CALC_NUM_ERR -eq 1 ] ; then
      # FIXME: Must remove commas here because I cannot figure out how to pass
      # the regular expression to ndiff.awk FS without the interpreter giving
      # this error for FS='[ \t,()]':
      # awk: fatal: Invalid regular expression: /'[/
      $CPPTRAJ_NCDUMP -n nctest $F1 | grep -v "==>\|:programVersion" | sed 's/,/ /g' > nc0.save
      $CPPTRAJ_NCDUMP -n nctest $F2 | grep -v "==>\|:programVersion" | sed 's/,/ /g' > nc0
    else
      $CPPTRAJ_NCDUMP -n nctest $F1 | grep -v "==>\|:programVersion" > nc0.save
      $CPPTRAJ_NCDUMP -n nctest $F2 | grep -v "==>\|:programVersion" > nc0
    fi
    DoTest $DIFFARGS 
    $CPPTRAJ_RM nc0.save nc0
  fi
}

# ------------------------------------------------------------------------------
# GetResultsFiles() <mode> <base>
#   Consolidate files with name <base>. If <mode> is 'single' assume one file in
#   the current directory, otherwise assume files are in subdirectories.
GetResultsFiles() {
  if [ ! -z "$TEST_DIRS" ] ; then
    RESULTSFILES=''
    for DIR in $TEST_DIRS ; do
      if [ -f "$DIR/$2" ] ; then
        RESULTSFILES="$RESULTSFILES $DIR/$2"
      fi
    done
  else
    if [ "$1" = 'single' ] ; then
      RESULTSFILES=`ls $2 2> /dev/null`
    else
      RESULTSFILES=`ls */$2 2> /dev/null`
    fi
  fi
}

# ------------------------------------------------------------------------------
ParseValgrindOut() {
  awk 'BEGIN{
    n_vg_err = 0;
    n_vg_allocs = 0;
    n_vg_frees = 0;
    n_vg_allocated = 0;
    n_vg_leaks = 0;
    currentTest = "";
    ntests = 0;
  }{
    if ($1 != currentTest) {
      ntests++;
      currentTest = $1;
    }
    if ($2 == "total" && $3 == "heap") {
      n_vg_allocs += $5;
      n_vg_frees += $7;
      n_vg_allocated += $9;
    } else if ($2 == "ERROR" && $3 == "SUMMARY:") {
      n_vg_err += $4;
    } else if ($2 == "LEAK" && $3 == "SUMMARY:") {
      n_vg_leaks++;
    }
  }END{
    printf("Valgrind summary for %i program executions:\n", ntests);
    printf("    %i errors.\n", n_vg_err);
    printf("    %i allocs, %i frees, %i bytes allocated.\n",
           n_vg_allocs, n_vg_frees, n_vg_allocated);
    printf("    %i memory leaks.\n", n_vg_leaks);
  }' $1
}

# ------------------------------------------------------------------------------
# Summary()
#  Print a summary of results in all CPPTRAJ_TEST_RESULTS files. Optionally
#  print an error summary and/or valgrind summary.
Summary() {
  MODE=''
  if [ "$CPPTRAJ_TEST_MODE" = 'master' ] ; then
    MODE='multi'
  else
    MODE='single'
  fi
  if [ ! -z "$CPPTRAJ_TEST_RESULTS" ] ; then
    GetResultsFiles $MODE $CPPTRAJ_TEST_RESULTS
    #echo "DEBUG: Getting results from $RESULTSFILES"
    if [ ! -z "$RESULTSFILES" ] ; then
      if [ "$MODE" = 'multi' ] ; then
        cat $RESULTSFILES > $CPPTRAJ_TEST_RESULTS
      fi
      echo "===================== TEST SUMMARY ======================"
      awk 'BEGIN{
        comparisons_ok = 0;   # Number of OK comparisons
        comparisons_warn = 0; # Number of comparisons with a warning
        comparisons_diff = 0; # Number of comparisons with diff/error
        program_err = 0;      # Number of program errors.
        program_exe = 0;      # Number of program executions.
        ntests = 0;           # Number of tests.
        ok_tests = 0;         # Number of completely OK tests.
      }{
        if ( $1 == "TEST:" )
          ntests++;
        else if ($1 == "Warning:")
          comparisons_warn++;
        else if ($1 == "CPPTRAJ:")
          program_exe++;
        else if ($1 == "Error:")
          program_err++;
        else if ($NF == "OK.")
          comparisons_ok++;
        else if ($NF == "different.")
          comparisons_diff++;
        else if ($(NF-1) == "not" && $NF == "found.")
          comparisons_diff++;
        else if ($(NF-1) == "comparisons" && $NF == "passed.")
          ok_tests++;
      }END{
        n_comps = comparisons_ok + comparisons_warn + comparisons_diff;
        printf("  %i out of %i comparisons OK (%i failed, %i warnings).\n",
               comparisons_ok, n_comps, comparisons_diff, comparisons_warn);
        printf("  %i out of %i program executions completed.\n",
               program_exe - program_err, program_exe);
        printf("  %i out of %i tests completed with no issues.\n",
               ok_tests, ntests);
      }' $CPPTRAJ_TEST_RESULTS 
    else
      echo "No test results files ($CPPTRAJ_TEST_RESULTS) found."
    fi
  fi
  # Error summary
  if [ ! -z "$CPPTRAJ_TEST_ERROR" ] ; then
    GetResultsFiles $MODE $CPPTRAJ_TEST_ERROR
    #echo "DEBUG: Getting errors from $RESULTSFILES"
    if [ ! -z "$RESULTSFILES" ] ; then
      if [ "$MODE" = 'multi' ] ; then
        cat $RESULTSFILES > $CPPTRAJ_TEST_ERROR
      fi
      if [ $SHOWERRORS -eq 1 ] ; then
        echo ""
        echo "Obtained the following errors:"
        echo "---------------------------------------------------------"
        cat $CPPTRAJ_TEST_ERROR
        echo "---------------------------------------------------------"
      fi
    fi
  fi
}

# ------------------------------------------------------------------------------
# RunCpptraj() <title>
#   Run cpptraj test with given title.
RunCpptraj() {
  # If only cleaning requested no run needed, exit now
  if [ $CPPTRAJ_TEST_CLEAN -eq 1 ] ; then
    exit 0
  fi
  ((PROGCOUNT++))
  echo ""
  echo "  CPPTRAJ: $1"
  if [ -z "$CPPTRAJ_DACDIF" ] ; then
    OUT "  CPPTRAJ: $1"
  fi
  if [ ! -z "$CPPTRAJ_DEBUG" ] ; then
    echo "$CPPTRAJ_TIME $DO_PARALLEL $VALGRIND $CPPTRAJ $TOP $INPUT $CPPTRAJ_DEBUG >> $CPPTRAJ_OUTPUT 2>>$CPPTRAJ_ERROR"
  fi
  $CPPTRAJ_TIME $DO_PARALLEL $VALGRIND $CPPTRAJ $TOP $INPUT $CPPTRAJ_DEBUG>> $CPPTRAJ_OUTPUT 2>>$CPPTRAJ_ERROR
  STATUS=$?
  #echo "DEBUG: Cpptraj exited with status $STATUS"
  if [ $STATUS -ne 0 ] ; then
    if [ -z "$CPPTRAJ_DACDIF" ] ; then
      echo "Error: cpptraj exited with status $STATUS" > /dev/stderr
      OutBoth "Error: cpptraj exited with status $STATUS"
    else
      echo "$CPPTRAJ: Program error"
    fi
    ((PROGERROR++))
  fi
}

# ------------------------------------------------------------------------------
# EndTest()
#   Print a summary of the current tests results. Should be called at the end of
#   every test script.
EndTest() {
  #echo "DEBUG: EndTest"
  # Report only when not using dacdif 
  if [ -z "$CPPTRAJ_DACDIF" ] ; then
    echo ""
    if [ $PROGERROR -gt 0 ] ; then
      echo "  $PROGERROR out of $PROGCOUNT executions exited with an error."
    fi
    if [ $ERRCOUNT -gt 0 ] ; then
      echo    "  $ERRCOUNT out of $NUMCOMPARISONS comparisons failed."
      OutBoth "  $ERRCOUNT out of $NUMCOMPARISONS comparisons failed."
    elif [ $WARNCOUNT -gt 0 ] ; then
      ((PASSCOUNT = $NUMCOMPARISONS - $WARNCOUNT))
      echo "  $PASSCOUNT out of $NUMCOMPARISONS passing comparisons. $WARNCOUNT warnings."
      OUT  "  $PASSCOUNT out of $NUMCOMPARISONS passing comparisons. $WARNCOUNT warnings."
    else 
      echo "All $NUMCOMPARISONS comparisons passed." 
      OUT  "All $NUMCOMPARISONS comparisons passed."
    fi
    echo ""
    if [ ! -z "$VALGRIND" ] ; then
      echo "Valgrind summary:"
      grep ERROR $CPPTRAJ_ERROR
      grep heap $CPPTRAJ_ERROR
      grep LEAK $CPPTRAJ_ERROR
      echo ""
      echo "Valgrind summary:" >> $CPPTRAJ_TEST_RESULTS
      grep ERROR $CPPTRAJ_ERROR >> $CPPTRAJ_TEST_RESULTS
      grep heap $CPPTRAJ_ERROR >> $CPPTRAJ_TEST_RESULTS
      grep LEAK $CPPTRAJ_ERROR >> $CPPTRAJ_TEST_RESULTS
    fi
    if [ $CPPTRAJ_PROFILE -eq 1 ] ; then
      if [ -e 'gmon.out' ] ; then
        gprof $CPPTRAJ > profiledata.txt
      fi
    fi
  fi
  if [ $SUMMARY -ne 0 ] ; then
    # NOTE: SUMMARY should only be set here if single test.
    Summary
  fi
}

# ------------------------------------------------------------------------------
# CleanFiles() <file1> ... <fileN>
#   For every arg passed to the function, check for the file or directory and
#   remove it.
CleanFiles() {
  while [ ! -z "$1" ] ; do
    if [ -d "$1" ] ; then
      rmdir $1
    elif [ -f "$1" ] ; then
      $CPPTRAJ_RM $1
    fi
    shift
  done
  # If only cleaning requested no run needed, exit now
  if [ $CPPTRAJ_TEST_CLEAN -eq 1 ] ; then
    exit 0
  fi
}

#-------------------------------------------------------------------------------
# NotParallel() <Test title>
#   Used to indicate a test should not be run in parallel.
NotParallel() {
  if [ ! -z "$DO_PARALLEL" ] ; then
    echo ""
    echo "  CPPTRAJ: $1"
    echo "  This test cannot be run in parallel. Skipping test."
    return 1
 fi
 return 0
}

#-------------------------------------------------------------------------------
# SetNthreads()
#   Use CPPTRAJ_NPROC to set N_THREADS if not already set.
SetNthreads() {
  if [ -z "$N_THREADS" ] ; then
    if [ ! -f "$CPPTRAJ_NPROC" ] ; then
      export N_THREADS=1
      return 1
    fi
    export N_THREADS=`$DO_PARALLEL $CPPTRAJ_NPROC`
    echo "  $N_THREADS MPI threads."
  fi
  return 0
}

#-------------------------------------------------------------------------------
# RequiresThreads() <# threads> <Test title>
RequiresThreads() {
  if [ ! -z "$DO_PARALLEL" ] ; then
    SetNthreads
    if [ $? -ne 0 ] ; then
      echo "Error: Program to find # threads not found ($CPPTRAJ_NPROC)" > /dev/stderr
      echo "Error: Test requires $1 parallel threads. Attempting to run test anyway." > /dev/stderr
      return 0
    fi
    REMAINDER=`echo "$N_THREADS % $1" | bc`
    if [ -z "$REMAINDER" -o $REMAINDER -ne 0 ] ; then
      echo ""
      if [ ! -z "$2" ] ; then
        echo "  CPPTRAJ: $2"
      fi
      echo "  Warning: Test requires a multiple of $1 parallel threads. Skipping."
      return 1
    fi
  fi
  return 0
}

#-------------------------------------------------------------------------------
# MaxThreads() <# threads> <Test title>
MaxThreads() {
  if [ ! -z "$DO_PARALLEL" ] ; then
    SetNthreads
    if [ $? -ne 0 ] ; then
      echo "Error: Program to find # threads not found ($CPPTRAJ_NPROC)" > /dev/stderr
      echo "Error: Test can only run with $1 or fewer threads. Attempting to run test anyway." > /dev/stderr
      return 0
    fi
    if [ $N_THREADS -gt $1 ] ; then
      echo ""
      if [ ! -z "$2" ] ; then
        echo "  CPPTRAJ: $2"
      fi
      echo "  Warning: Test can only run with $1 or fewer parallel threads. Skipping."
      return 1
    fi
  fi
  return 0
}

#-------------------------------------------------------------------------------
# Help(): Print help
Help() {
  echo "Command line flags"
  if [ "$CPPTRAJ_TEST_MODE" = 'master' ] ; then
    echo "  <test dir 1> [<test dir 2> ...]"
  fi
  echo "  mpi             : Use MPI version of CPPTRAJ (automatically used if DO_PARALLEL set)."
  echo "  openmp          : Use OpenMP version of CPPTRAJ."
  echo "  cuda            : Use CUDA version of CPPTRAJ."
  echo "  time            : Time the test with 'time'."
  echo "  stdout          : Print CPPTRAJ test output to STDOUT."
  echo "  vg              : Run test with 'valgrind' memcheck."
  echo "  vgh             : Run test with 'valgrind' helgrind."
  echo "  clean           : Do not run test; clean test output."
  echo "  summary         : Print summary of tests run."
  echo "  showerrors      : (summary only) Print all test errors to STDOUT after summary."
  echo "  -nodacdif       : Do not use dacdif for test comparisons (AmberTools only)."
  echo "  -d              : Run CPPTRAJ with global debug level 4."
  echo "  -debug <#>      : Run CPPTRAJ with global debug level #."
  echo "  -cpptraj <file> : Use CPPTRAJ binary <file>."
  echo "  -ambpdb <file>  : Use AMBPDB binary <file>."
  echo "  -profile        : Profile results with 'gprof' (requires special compile)."
  echo "Important environment variables"
  echo "  DO_PARALLEL     : MPI run command."
  echo "  N_THREADS       : Number of MPI threads (only needed if 'DO_PARALLEL nproc' fails)."
  echo "  OMP_NUM_THREADS : Number of OpenMP threads to use."
  echo "  CPPTRAJ_DEBUG   : Can be used to pass extra flags to CPPTRAJ."
  echo ""
}

#-------------------------------------------------------------------------------
# CmdLineOpts()
#   Process test script command line options. Only executed if
#   CPPTRAJ_TEST_SETUP is not already set.
CmdLineOpts() {
  CPPTRAJ_TEST_CLEAN=0 # Will be exported
  CPPTRAJ_PROFILE=0    # Will be exported
  SFX_OMP=0
  SFX_CUDA=0
  SFX_MPI=0
  GET_TIMING=0
  while [ ! -z "$1" ] ; do
    case "$1" in
      "clean"     ) CPPTRAJ_TEST_CLEAN=1 ;;
      "summary"   ) SUMMARY=1 ;;
      "showerrors") SHOWERRORS=1 ;;
      "stdout"    ) CPPTRAJ_OUTPUT='/dev/stdout' ;;
      "openmp"    ) SFX_OMP=1 ;;
      "cuda"      ) SFX_CUDA=1 ;;
      "mpi"       ) SFX_MPI=1 ;;
      "vg"        ) VGMODE=1 ;;
      "vgh"       ) VGMODE=2 ;;
      "time"      ) GET_TIMING=1 ;;
      "-d"        ) CPPTRAJ_DEBUG="$CPPTRAJ_DEBUG -debug 4" ;;
      "-debug"    ) shift ; CPPTRAJ_DEBUG="$CPPTRAJ_DEBUG -debug $1" ;;
      "-nodacdif" ) USE_DACDIF=0 ;;
      "-cpptraj"  ) shift ; export CPPTRAJ=$1 ; echo "Using cpptraj: $CPPTRAJ" ;;
      "-ambpdb"   ) shift ; export AMBPDB=$1  ; echo "Using ambpdb: $AMBPDB" ;;
      "--target"  ) shift ; TARGET=$1 ;;
     "-profile"   ) CPPTRAJ_PROFILE=1 ; echo "Performing gnu profiling during EndTest." ;;
      "-h" | "--help" ) Help ; exit 0 ;;
      "alltests"  )
        echo "Running all tests in Test_* directories."
        for DIR in `ls -d Test_*` ; do
          if [ -f "$DIR/RunTest.sh" ] ; then
            TEST_DIRS="$TEST_DIRS $DIR"
          fi
        done
        ;;
      *           )
        if [ -d "$1" -a -f "$1/RunTest.sh" ] ; then
          # Assume this is a test we want to run.
          TEST_DIRS="$TEST_DIRS $1"
        else
          echo "Error: Unknown option: $1" > /dev/stderr
          exit 1
        fi
        ;;
    esac
    shift
  done
  if [ $CPPTRAJ_PROFILE -eq 1 ] ; then
    Required "gprof"
  fi
  export CPPTRAJ_TEST_CLEAN
  export CPPTRAJ_DEBUG
  export CPPTRAJ_PROFILE
  # If DO_PARALLEL has been set force MPI
  if [ ! -z "$DO_PARALLEL" ] ; then
    SFX_MPI=1
  elif [ $SFX_MPI -eq 1 ] ; then
    echo "Error: 'mpi' specified but DO_PARALLEL not set." > /dev/stderr
    exit 1
  fi
  # Warn if using OpenMP but OMP_NUM_THREADS not set.
  if [ $SFX_OMP -eq 1 -a -z "$OMP_NUM_THREADS" ] ; then
    echo "Warning: Using OpenMP but OMP_NUM_THREADS is not set."
  fi
  # Set up SFX
  if [ "$SFX_MPI"  -eq 1 ] ; then SFX="$SFX.MPI"  ; fi
  if [ "$SFX_OMP"  -eq 1 ] ; then SFX="$SFX.OMP"  ; fi
  if [ "$SFX_CUDA" -eq 1 ] ; then SFX="$SFX.cuda" ; fi
}

#-------------------------------------------------------------------------------
# Required() <binary>
#   Insure that specified binary is in the PATH
Required() {
  if [ -z "`which $1`" ] ; then
    echo "Error: Required binary '$1' not found." > /dev/stderr
    exit 1
  fi
}

#-------------------------------------------------------------------------------
# CheckDefines()
#   Check how CPPTRAJ was compiled.
CheckDefines() {
  for DEFINE in `$CPPTRAJ --defines` ; do
    case "$DEFINE" in
      '-DHASGZ'         ) export CPPTRAJ_ZLIB=$DEFINE ;;
      '-DHASBZ2'        ) export CPPTRAJ_BZLIB=$DEFINE ;;
      '-DBINTRAJ'       ) export CPPTRAJ_NETCDFLIB=$DEFINE ;;
      '-DMPI'           ) export CPPTRAJ_MPILIB=$DEFINE ;;
      '-DNO_MATHLIB'    ) export CPPTRAJ_NOMATHLIB=$DEFINE ;;
      '-D_OPENMP'       ) export CPPTRAJ_OPENMP=$DEFINE ;;
      '-DHAS_PNETCDF'   ) export CPPTRAJ_PNETCDFLIB=$DEFINE ;;
      '-DUSE_SANDERLIB' ) export CPPTRAJ_SANDERLIB=$DEFINE ;;
      '-DCUDA'          ) export CPPTRAJ_CUDA=$DEFINE ;;
      '-DNO_XDRFILE'    ) export CPPTRAJ_NO_XDRFILE=$DEFINE ;;
    esac
  done
  #echo "DEBUG: $ZLIB $BZLIB $NETCDFLIB $MPILIB $NOMATHLIB $OPENMP $PNETCDFLIB $SANDERLIB $CUDA $NO_XDRFILE"
}

#-------------------------------------------------------------------------------
# SetBinaries()
#   Set paths for all binaries if not already set.
SetBinaries() {
  # Guess where things might be depending on if we are in AT or not, etc.
  if [ $STANDALONE -eq 0 ] ; then
    DIRPREFIX=$AMBERHOME
  elif [ -z "$CPPTRAJHOME" ] ; then
    DIRPREFIX=../../
  else
    DIRPREFIX=$CPPTRAJHOME
  fi
  # Determine location of CPPTRAJ if not specified.
  if [ -z "$CPPTRAJ" ] ; then
    CPPTRAJ=$DIRPREFIX/bin/cpptraj$SFX
    if [ ! -f "$CPPTRAJ" ] ; then
      echo "Error: CPPTRAJ $CPPTRAJ not found." > /dev/stderr
      exit 1
    fi
    export CPPTRAJ
  fi
  # Check how cpptraj was compiled.
  CheckDefines
  # If compiled with NetCDF support ensure ncdump is available
  if [ ! -z "$CPPTRAJ_NETCDFLIB" ] ; then
    Required "ncdump"
    export CPPTRAJ_NCDUMP=`which ncdump`
  fi
  # Set up valgrind if necessary.
  if [ $VGMODE -ne 0 ] ; then
    echo "  Using valgrind."
    CPPTRAJ_ERROR='valgrind.out'
    if [ ! -z "$CPPTRAJ_CUDA" ] ; then
      if [ $VGMODE -eq 1 ] ; then
        Required "cuda-memcheck"
        VALGRIND='cuda-memcheck'
      else
        echo "Error: 'helgrind' not supported for CUDA." > /dev/stderr
        exit 1
      fi
    else
      Required "valgrind"
      if [ $VGMODE -eq 1 ] ; then
        VALGRIND="valgrind --tool=memcheck --leak-check=yes --show-reachable=yes"
      elif [ $VGMODE -eq 2 ] ; then
        VALGRIND="valgrind --tool=helgrind"
      else
        echo "Error: Unsupported VGMODE $VGMODE" > /dev/stderr
        exit 1
      fi
    fi
    export VALGRIND
  fi
  # Determine location of AMBPDB if not specified.
  if [ -z "$AMBPDB" ] ; then
    AMBPDB=$DIRPREFIX/bin/ambpdb
    if [ ! -f "$AMBPDB" ] ; then
      echo "Warning: AMBPDB not present."
    fi
    export AMBPDB
  fi
  # Determine location of ndiff.awk
  if [ -z "$CPPTRAJ_NDIFF" ] ; then
    if [ $STANDALONE -eq 0 ] ; then
      CPPTRAJ_NDIFF=$DIRPREFIX/test/ndiff.awk
    else
      CPPTRAJ_NDIFF=$DIRPREFIX/util/ndiff/ndiff.awk
    fi
    if [ ! -f "$CPPTRAJ_NDIFF" ] ; then
      echo "Error: 'ndiff.awk' not present: $CPPTRAJ_NDIFF" > /dev/stderr
      exit 1
    fi
    export CPPTRAJ_NDIFF
  fi
  # Determine location of nproc/numprocs
  if [ -z "$CPPTRAJ_NPROC" ] ; then
    if [ ! -z "$DO_PARALLEL" ] ; then
      if [ $STANDALONE -eq 0 ] ; then
        CPPTRAJ_NPROC=$DIRPREFIX/AmberTools/test/numprocs
      else
        CPPTRAJ_NPROC=$DIRPREFIX/test/nproc
      fi
      if [ -z "$CPPTRAJ_NPROC" ] ; then
        echo "Error: nproc $CPPTRAJ_NPROC not found." > /dev/stderr
        exit 1
      fi
      export CPPTRAJ_NPROC
      SetNthreads
    fi
  fi
  # Determine timing if necessary
  if [ $GET_TIMING -eq 1 ] ; then
    Required "time"
    export CPPTRAJ_TIME=`which time`
  fi
  # Report binary details
  if [ $STANDALONE -eq 1 ] ; then
    ls -l $CPPTRAJ
    if [ ! -z "$AMBPDB" ] ; then
      ls -l $AMBPDB
    fi
  fi
  # Print DEBUG info
  if [ ! -z "$CPPTRAJ_DEBUG" ] ; then
    if [ $STANDALONE -eq 1 ] ; then
      echo "DEBUG: Standalone mode."
    else
      echo "DEBUG: AmberTools mode."
    fi
    echo "DEBUG: CPPTRAJ: $CPPTRAJ"
    echo "DEBUG: AMBPDB:  $AMBPDB"
    echo "DEBUG: NPROC:   $CPPTRAJ_NPROC"
    echo "DEBUG: NCDUMP:  $CPPTRAJ_NCDUMP"
    echo "DEBUG: DIFFCMD: $CPPTRAJ_DIFF"
    echo "DEBUG: DACDIF:  $CPPTRAJ_DACDIF"
    echo "DEBUG: NDIFF:   $CPPTRAJ_NDIFF"
  fi
}

# ------------------------------------------------------------------------------
# Library Checks - Tests that depend on certain libraries like Zlib can run
# these to make sure cpptraj was compiled with that library - exit gracefully
# if not.
# Should not be called if CLEAN==1, CleanFiles should always be called first.
CheckZlib() {
  if [ -z "$CPPTRAJ_ZLIB" ] ; then
    echo "This test requires zlib. Cpptraj was compiled without zlib support."
    echo "Skipping test."
    exit 0
  fi
}

CheckBzlib() {
  if [ -z "$CPPTRAJ_BZLIB" ] ; then
    echo "This test requires bzlib. Cpptraj was compiled without bzlib support."
    echo "Skipping test."
    exit 0
  fi
}

CheckNetcdf() {
  if [ -z "$CPPTRAJ_NETCDFLIB" ] ; then
    echo "This test requires NetCDF. Cpptraj was compiled without NetCDF support."
    echo "Skipping test."
    exit 0
  fi
}

CheckPtrajAnalyze() {
  if [ ! -z "$CPPTRAJ_NOMATHLIB" ] ; then
    echo "This test requires LAPACK/ARPACK/BLAS routines."
    echo "Cpptraj was compiled with -DNO_MATHLIB. Skipping test."
    exit 0
  fi
}

CheckSanderlib() {
  if [ -z "$CPPTRAJ_SANDERLIB" ] ; then
    if [ ! -z "$1" ] ; then
      DESCRIP=$1
    else
      DESCRIP="This test"
    fi
    echo "$DESCRIP requires compilation with the Sander API from AmberTools."
    echo "Skipping test."
    return 1
  fi
  return 0
}

CheckPnetcdf() {
  if [ ! -z "$DO_PARALLEL" ] ; then
    DESCRIP="This test"
    if [ ! -z "$1" ] ; then
      DESCRIP="Test '$1'"
    fi
    if [ -z "$CPPTRAJ_PNETCDFLIB" ] ; then
      echo "$DESCRIP requires compilation with parallel NetCDF."
      echo "Cpptraj was compiled without parallel NetCDF support. Skipping test."
      return 1
    fi
  fi
  return 0
}

# ==============================================================================
# M A I N
# ==============================================================================
#echo "DEBUG: Begin MasterTest.sh. $*"
#echo "DEBUG: CPPTRAJ_TEST_MODE: $CPPTRAJ_TEST_MODE"
if [ -z "$CPPTRAJ_TEST_SETUP" ] ; then
  echo "DEBUG: Initial test setup."
  # MasterTest.sh has not been called yet; set up test environment.
  export CPPTRAJ_TEST_ROOT=`pwd`
  # Ensure required binaries are set up
  if [ -z "$CPPTRAJ_RM" ] ; then
    # TODO is this being too paranoid?
    if [ ! -f '/bin/rm' ] ; then
      echo "Error: Required binary '/bin/rm' not found." > /dev/stderr
      exit 1
    fi
    export CPPTRAJ_RM='/bin/rm -f'
  fi
  Required "grep"
  Required "sed"
  Required "awk"
  Required "rmdir"
  # Set some defaults
  export CPPTRAJ_TEST_RESULTS='Test_Results.dat'
  export CPPTRAJ_TEST_ERROR='Test_Error.dat'
  CPPTRAJ_OUTPUT='test.out'
  CPPTRAJ_ERROR='/dev/stderr'
  # Process command line options
  CmdLineOpts $*
  # If not cleaning see what else needs to be set up. 
  if [ $CPPTRAJ_TEST_CLEAN -eq 0 ] ; then
    # Determine standalone or AmberTools
    if [ ! -z "`echo "$CPPTRAJ_TEST_ROOT" | grep AmberTools`" ] ; then
      # Assume AmberTools. Need AMBERHOME.
      STANDALONE=0
      if [ -z "$AMBERHOME" ] ; then
        echo "Error: In AmberTools and AMBERHOME is not set. Required for tests." > /dev/stderr
        exit 1
      fi
    else
      # Standalone. Never use dacdif.
      STANDALONE=1
      USE_DACDIF=0
    fi
    # Determine if diff or dacdif will be used.
    CPPTRAJ_DIFF=''
    CPPTRAJ_DACDIF=''
    if [ $USE_DACDIF -eq 1 ] ; then
      CPPTRAJ_DACDIF="$AMBERHOME/test/dacdif"
      if [ ! -f "$CPPTRAJ_DACDIF" ] ; then
        echo "$CPPTRAJ_DACDIF not found. Required for AmberTools tests."
        exit 1
      fi
    else
      Required "diff"
      CPPTRAJ_DIFF=`which diff`
    fi
    export CPPTRAJ_DACDIF
    export CPPTRAJ_DIFF
    # Determine binary locations
    SetBinaries
    # If CPPTRAJ_TEST_OS is not set, assume linux. FIXME needed?
    if [ -z "$CPPTRAJ_TEST_OS" ] ; then
      export CPPTRAJ_TEST_OS='linux'
    fi
    if [ ! -z "$DIFFOPTS" ] ; then
      echo "Warning: DIFFOPTS is set to '$DIFFOPTS'"
    fi
  fi # END if not cleaning
  # Export test output and error file names
  export CPPTRAJ_OUTPUT
  export CPPTRAJ_ERROR
  # Initial setup complete
  #echo "DEBUG: Initial test setup complete."
  export CPPTRAJ_TEST_SETUP='yes'
fi # END initial setup

# Determine mode of execution: individual test or multiple tests.
if [ "$CPPTRAJ_TEST_MODE" = 'master' ] ; then
  # Executed from CpptrajTest.sh. Assume we are executing multiple
  # tests.
  #echo "DEBUG: Executing multiple tests."
  if [ ! -z "$TEST_DIRS" ] ; then
    #echo "DEBUG: Running tests in specified directories."
    for DIR in $TEST_DIRS ; do
      cd $CPPTRAJ_TEST_ROOT/$DIR && ./RunTest.sh
    done
    if [ $SUMMARY -ne 0 ] ; then
      cd $CPPTRAJ_TEST_ROOT
      Summary
    fi
  else
    # Probably executed via make.
    # If summary requested just do that and exit.
     if [ $SUMMARY -ne 0 ] ; then
      Summary
      exit 0
    fi
    # Need a Makefile.
    if [ ! -f 'Makefile' ] ; then
      echo "Error: test Makefile not found." > /dev/stderr
      exit 1
    fi
    Required "make"
    if [ -z "$TARGET" ] ; then
      # If no target specified, probably not executed via make.
      echo "Error: No test directories specified." > /dev/stderr
      exit 1
    fi
    make $TARGET
    if [ $? -ne 0 ] ; then
      exit 1
    fi
  fi
else
  # Assume we are only executing a single test.
  #echo "DEBUG: Executing single test."
  # Single test.
  # Always clean up individual test output and error files
  if [ "$CPPTRAJ_OUTPUT" != '/dev/stdout' -a -f "$CPPTRAJ_OUTPUT" ] ; then
    $CPPTRAJ_RM $CPPTRAJ_OUTPUT
  fi
  if [ "$CPPTRAJ_ERROR" != '/dev/stderr' -a -f "$CPPTRAJ_ERROR" ] ; then
    $CPPTRAJ_RM $CPPTRAJ_ERROR
  fi
  if [ -f "$CPPTRAJ_TEST_RESULTS" ] ; then
    $CPPTRAJ_RM $CPPTRAJ_TEST_RESULTS
  fi
  if [ -f "$CPPTRAJ_TEST_ERROR" ] ; then
    $CPPTRAJ_RM $CPPTRAJ_TEST_ERROR
  fi
  if [ -f 'valgrind.out' ] ; then
    $CPPTRAJ_RM valgrind.out
  fi
  if [ -f 'test.out' ] ; then
    $CPPTRAJ_RM test.out
  fi
  THREADFILES=`ls Thread.* 2> /dev/null`
  if [ ! -z "$THREADFILES" ] ; then
    for FILE in $THREADFILES ; do
      $CPPTRAJ_RM $FILE
    done
  fi
  if [ "$CPPTRAJ_TEST_CLEAN" -eq 0 ] ; then
    TEST_WORKDIR=`pwd`
    TestHeader '/dev/stdout'
    if [ -z "$CPPTRAJ_DACDIF" ] ; then
      TestHeader "$CPPTRAJ_TEST_RESULTS"
    fi
  fi
fi

#echo "DEBUG: End MasterTest.sh."
