description="""
An app to calculate the correction for releasing Cartesian restraints
"""
from Sire.Tools import CartesianCorrection
from Sire.Tools import readParams

import argparse
import sys
import os

parser = argparse.ArgumentParser(description="""Calculates the correction for releasing Cartesian restraints,
                                             accounting for the standard state.""",
                                epilog="cartesian_correction is built using Sire and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org",
                                 prog="cartesian_correction")

parser.add_argument('-C', '--config', nargs="?",
                    help='A config file used for the simulation (all are suitable)')

parser.add_argument('-H', '--help-config', action="store_true",
                    help="Get additional help regarding all of the parameters "
                         "(and their default values) that can be "
                         "set in the optionally-supplied CONFIG file")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")
                    
parser.add_argument('--verbose', action="store_true",
                    help="Print debug output.")

sys.stdout.write("\n")
args = parser.parse_args()

must_exit = False

if args.author:
    print("""\ncartesian_correction was written by Finlay Clark (C) 2022, but is very closely
     based on standardstatecorrection, written by Julien Michel and Stefano Bosisio (C) 2017""")
    must_exit = True

if args.version:
    print("cartesian_correction -- from Sire release version <%s>" %Sire.__version__)
    print("This particular release can be downloaded here: "
          "https://github.com/michellab/Sire/releases/tag/v%s" %Sire.__version__)
    must_exit = True

if args.help_config:
    OpenMMMD.Parameter.printAll(True)
    must_exit = True

if must_exit:
    sys.exit(0)

if args.config:
    config_file = args.config
    if not os.path.exists(config_file):
        print("(cannot find configuration file %s)" % config_file)
        sys.exit(-1)
    else:
        print("Loading temperature and restraints information from file %s" % args.config)
        params = readParams(args.config)
else:
    print("Please supply a configuration file")
    sys.exit(-1)

if args.verbose:
    params["verbose"] = True

print("\nCalculating correction for Cartesian restraints using temperature and"
      "restraint information from %s." % args.config)

CartesianCorrection.run(params)


