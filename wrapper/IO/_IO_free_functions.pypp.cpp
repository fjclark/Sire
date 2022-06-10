// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "_IO_free_functions.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/getinstalldir.h"

#include "SireError/errors.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/mgname.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidx.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "biosimspace.h"

#include "moleculeparser.h"

#include "biosimspace.h"

#include "SireBase/getinstalldir.h"

#include "SireError/errors.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/mgname.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidx.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "biosimspace.h"

#include "moleculeparser.h"

#include "biosimspace.h"

#include "SireBase/getinstalldir.h"

#include "SireError/errors.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/mgname.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidx.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "biosimspace.h"

#include "moleculeparser.h"

#include "biosimspace.h"

#include "SireBase/getinstalldir.h"

#include "SireError/errors.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/mgname.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidx.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "biosimspace.h"

#include "moleculeparser.h"

#include "biosimspace.h"

#include "SireBase/getinstalldir.h"

#include "SireError/errors.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/mgname.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidx.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "biosimspace.h"

#include "moleculeparser.h"

#include "biosimspace.h"

#include "SireBase/getinstalldir.h"

#include "SireError/errors.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/mgname.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidx.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "biosimspace.h"

#include "moleculeparser.h"

#include "biosimspace.h"

#include "SireBase/getinstalldir.h"

#include "SireError/errors.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/mgname.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidx.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "biosimspace.h"

#include "moleculeparser.h"

#include "biosimspace.h"

#include "SireBase/getinstalldir.h"

#include "SireError/errors.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/mgname.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidx.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "biosimspace.h"

#include "moleculeparser.h"

#include "biosimspace.h"

#include "SireBase/getinstalldir.h"

#include "SireError/errors.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/mgname.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidx.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "biosimspace.h"

#include "moleculeparser.h"

#include "biosimspace.h"

#include "SireBase/getinstalldir.h"

#include "SireError/errors.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/mgname.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidx.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "biosimspace.h"

#include "moleculeparser.h"

#include "biosimspace.h"

#include "SireBase/getinstalldir.h"

#include "SireError/errors.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/mgname.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidx.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "biosimspace.h"

#include "moleculeparser.h"

#include "biosimspace.h"

#include "SireBase/getinstalldir.h"

#include "SireError/errors.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/mgname.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidx.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "biosimspace.h"

#include "moleculeparser.h"

#include "biosimspace.h"

#include "SireBase/getinstalldir.h"

#include "SireError/errors.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/connectivity.h"

#include "SireMol/mgname.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidx.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "biosimspace.h"

#include "moleculeparser.h"

#include "biosimspace.h"

void register_free_functions(){

    { //::SireIO::isAmberWater
    
        typedef bool ( *isAmberWater_function_type )( ::SireMol::Molecule const &,::SireBase::PropertyMap const & );
        isAmberWater_function_type isAmberWater_function_value( &::SireIO::isAmberWater );
        
        bp::def( 
            "isAmberWater"
            , isAmberWater_function_value
            , ( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() )
            , "Test whether the passed water molecule matches standard AMBER\nformat water topologies.\n\nPar:am molecule\nThe molecule to test.\n\nPar:am map\nA dictionary of user-defined molecular property names.\n\nRetval: is_water\nWhether the molecule is an AMBER format water.\n" );
    
    }

    { //::SireIO::isGromacsWater
    
        typedef bool ( *isGromacsWater_function_type )( ::SireMol::Molecule const &,::SireBase::PropertyMap const & );
        isGromacsWater_function_type isGromacsWater_function_value( &::SireIO::isGromacsWater );
        
        bp::def( 
            "isGromacsWater"
            , isGromacsWater_function_value
            , ( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() )
            , "Test whether the passed water molecule matches standard GROMACS\nformat water topologies.\n\nPar:am molecule\nThe molecule to test.\n\nPar:am map\nA dictionary of user-defined molecular property names.\n\nRetval: is_water\nWhether the molecule is a GROMACS format water.\n" );
    
    }

    { //::SireIO::isWater
    
        typedef bool ( *isWater_function_type )( ::SireMol::Molecule const &,::SireBase::PropertyMap const & );
        isWater_function_type isWater_function_value( &::SireIO::isWater );
        
        bp::def( 
            "isWater"
            , isWater_function_value
            , ( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() )
            , "Test whether the passed water molecule matches standard water\ntopologies.\n\nPar:am molecule\nThe molecule to test.\n\nPar:am map\nA dictionary of user-defined molecular property names.\n\nRetval: is_water\nWhether the molecule is a water.\n" );
    
    }

    { //::SireIO::renumberConstituents
    
        typedef ::SireSystem::System ( *renumberConstituents_function_type )( ::SireSystem::System const &,unsigned int );
        renumberConstituents_function_type renumberConstituents_function_value( &::SireIO::renumberConstituents );
        
        bp::def( 
            "renumberConstituents"
            , renumberConstituents_function_value
            , ( bp::arg("system"), bp::arg("mol_offset")=(unsigned int)(0) )
            , "Renumber the constituents of a system (residues and atoms) so that\nthey are unique and are in ascending order.\n\nPar:am system\nThe molecular system of interest.\n\nPar:am mol_offset\nThe index of the molecule at which to begin renumbering.\n\nRetval: system\nThe system with renumbered constituents.\n" );
    
    }

    { //::SireIO::repartitionHydrogenMass
    
        typedef ::SireSystem::System ( *repartitionHydrogenMass_function_type )( ::SireSystem::System const &,double const,unsigned int const,::SireBase::PropertyMap const & );
        repartitionHydrogenMass_function_type repartitionHydrogenMass_function_value( &::SireIO::repartitionHydrogenMass );
        
        bp::def( 
            "repartitionHydrogenMass"
            , repartitionHydrogenMass_function_value
            , ( bp::arg("system"), bp::arg("factor")=4, bp::arg("water")=(unsigned int const)(0), bp::arg("map")=SireBase::PropertyMap() )
            , "Redistribute mass of heavy atoms connected to bonded hydrogens into\nthe hydrogen atoms. This allows use of larger simulation integration\ntime steps without encountering instabilities related to high-frequency\nhydrogen motion.\n\nPar:am system\nThe molecular system of interest.\n\nPar:am factor\nThe repartitioning scale factor. Hydrogen masses are scaled by\nthis amount.\n\nPar:am water\nWhether to repartiotion masses for water molecules:\n0 = yes, 1 = no, 2 = only water molecules.\n\nPar:am map\nA dictionary of user-defined molecular property names.\n\nRetval: system\nThe system with repartitioned hydrogen mass.\n" );
    
    }

    { //::SireIO::repartitionHydrogenMass
    
        typedef ::SireMol::Molecule ( *repartitionHydrogenMass_function_type )( ::SireMol::Molecule &,double const,unsigned int const,::SireBase::PropertyMap const & );
        repartitionHydrogenMass_function_type repartitionHydrogenMass_function_value( &::SireIO::repartitionHydrogenMass );
        
        bp::def( 
            "repartitionHydrogenMass"
            , repartitionHydrogenMass_function_value
            , ( bp::arg("molecule"), bp::arg("factor")=4, bp::arg("water")=(unsigned int const)(0), bp::arg("map")=SireBase::PropertyMap() )
            , "Redistribute mass of heavy atoms connected to bonded hydrogens into\nthe hydrogen atoms. This allows use of larger simulation integration\ntime steps without encountering instabilities related to high-frequency\nhydrogen motion.\n\nPar:am molecule\nThe molecule of interest.\n\nPar:am factor\nThe repartitioning scale factor. Hydrogen masses are scaled by\nthis amount.\n\nPar:am water\nWhether to repartiotion masses for water molecules:\n0 = yes, 1 = no, 2 = only water molecules.\n\nPar:am map\nA dictionary of user-defined molecular property names.\n\nRetval: system\nThe system with repartitioned hydrogen mass.\n" );
    
    }

    { //::SireIO::setAmberWater
    
        typedef ::SireSystem::System ( *setAmberWater_function_type )( ::SireSystem::System const &,::QString const &,::SireBase::PropertyMap const & );
        setAmberWater_function_type setAmberWater_function_value( &::SireIO::setAmberWater );
        
        bp::def( 
            "setAmberWater"
            , setAmberWater_function_value
            , ( bp::arg("system"), bp::arg("model"), bp::arg("map")=SireBase::PropertyMap() )
            , "Set all water molecules in the passed system to the appropriate AMBER\nformat topology.\n\nPar:am system\nThe molecular system of interest.\n\nPar:am model\nThe name of the water model.\n\nPar:am map\nA dictionary of user-defined molecular property names.\n\nRetval: system\nThe system with updated water topology.\n" );
    
    }

    { //::SireIO::setAmberWater
    
        typedef ::SireMol::SelectResult ( *setAmberWater_function_type )( ::SireMol::SelectResult const &,::QString const &,::SireBase::PropertyMap const & );
        setAmberWater_function_type setAmberWater_function_value( &::SireIO::setAmberWater );
        
        bp::def( 
            "setAmberWater"
            , setAmberWater_function_value
            , ( bp::arg("molecules"), bp::arg("model"), bp::arg("map")=SireBase::PropertyMap() )
            , "Set all water molecules in the passed system to the appropriate AMBER\nformat topology.\n\nPar:am system\nThe molecular system of interest.\n\nPar:am model\nThe name of the water model.\n\nPar:am map\nA dictionary of user-defined molecular property names.\n\nRetval: system\nThe system with updated water topology.\n" );
    
    }

    { //::SireIO::setGromacsWater
    
        typedef ::SireSystem::System ( *setGromacsWater_function_type )( ::SireSystem::System const &,::QString const &,::SireBase::PropertyMap const & );
        setGromacsWater_function_type setGromacsWater_function_value( &::SireIO::setGromacsWater );
        
        bp::def( 
            "setGromacsWater"
            , setGromacsWater_function_value
            , ( bp::arg("system"), bp::arg("model"), bp::arg("map")=SireBase::PropertyMap() )
            , "Set all water molecules in the passed system to the appropriate GROMACS\nformat topology.\n\nPar:am system\nThe molecular system of interest.\n\nPar:am model\nThe name of the water model.\n\nPar:am map\nA dictionary of user-defined molecular property names.\n\nRetval: system\nThe system with updated water topology.\n" );
    
    }

    { //::SireIO::setGromacsWater
    
        typedef ::SireMol::SelectResult ( *setGromacsWater_function_type )( ::SireMol::SelectResult const &,::QString const &,::SireBase::PropertyMap const & );
        setGromacsWater_function_type setGromacsWater_function_value( &::SireIO::setGromacsWater );
        
        bp::def( 
            "setGromacsWater"
            , setGromacsWater_function_value
            , ( bp::arg("molecules"), bp::arg("model"), bp::arg("map")=SireBase::PropertyMap() )
            , "Set all water molecules in the passed system to the appropriate GROMACS\nformat topology.\n\nPar:am system\nThe molecular system of interest.\n\nPar:am model\nThe name of the water model.\n\nPar:am map\nA dictionary of user-defined molecular property names.\n\nRetval: system\nThe system with updated water topology.\n\n" );
    
    }

    { //::SireIO::updateAndPreserveOrder
    
        typedef ::SireSystem::System ( *updateAndPreserveOrder_function_type )( ::SireSystem::System const &,::SireMol::Molecule const &,unsigned int );
        updateAndPreserveOrder_function_type updateAndPreserveOrder_function_value( &::SireIO::updateAndPreserveOrder );
        
        bp::def( 
            "updateAndPreserveOrder"
            , updateAndPreserveOrder_function_value
            , ( bp::arg("system"), bp::arg("molecule"), bp::arg("index") )
            , "Update a molecule in the system with a different UUID while\npreserving the molecular ordering. Normally we would need to\ndelete and re-add the molecule, which would place it at the\nend, even if the MolNum was unchanged.\n\nPar:am system\nThe molecular system of interest.\n\nPar:am molecule\nThe updated molecule.\n\nPar:am index\nThe index of the molecule in the system.\n\nRetval: system\nThe system with renumbered constituents.\n" );
    
    }

    { //::SireIO::updateCoordinatesAndVelocities
    
        typedef ::boost::tuples::tuple< SireSystem::System, QHash< SireMol::MolIdx, SireMol::MolIdx >, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( *updateCoordinatesAndVelocities_function_type )( ::SireSystem::System const &,::SireSystem::System const &,::QHash< SireMol::MolIdx, SireMol::MolIdx > const &,bool const,::SireBase::PropertyMap const &,::SireBase::PropertyMap const & );
        updateCoordinatesAndVelocities_function_type updateCoordinatesAndVelocities_function_value( &::SireIO::updateCoordinatesAndVelocities );
        
        bp::def( 
            "updateCoordinatesAndVelocities"
            , updateCoordinatesAndVelocities_function_value
            , ( bp::arg("system0"), bp::arg("system1"), bp::arg("molecule_mapping"), bp::arg("is_lambda1")=(bool const)(false), bp::arg("map0")=SireBase::PropertyMap(), bp::arg("map1")=SireBase::PropertyMap() )
            , "Update the coordinates and velocities of system0 with those from\nsystem1.\nPar:am system0\nThe reference system.\nPar:am system1\nThe updated system, where molecules may not be in the same order.\nPar:am map0\nA dictionary of user-defined molecular property names for system0.\nPar:am map1\nA dictionary of user-defined molecular property names for system1.\nRetval: system, mapping\nThe system with updated coordinates and velocities and a mapping\nbetween the molecule indices in both systems.\n" );
    
    }

    { //::SireIO::updateCoordinatesAndVelocities
    
        typedef ::boost::tuples::tuple< SireSystem::System, QHash< SireMol::MolIdx, SireMol::MolIdx >, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( *updateCoordinatesAndVelocities_function_type )( ::SireSystem::System const &,::SireSystem::System const &,::SireSystem::System const &,::QHash< SireMol::MolIdx, SireMol::MolIdx > const &,bool const,::SireBase::PropertyMap const &,::SireBase::PropertyMap const & );
        updateCoordinatesAndVelocities_function_type updateCoordinatesAndVelocities_function_value( &::SireIO::updateCoordinatesAndVelocities );
        
        bp::def( 
            "updateCoordinatesAndVelocities"
            , updateCoordinatesAndVelocities_function_value
            , ( bp::arg("original_system"), bp::arg("renumbered_system"), bp::arg("updated_system"), bp::arg("molecule_mapping"), bp::arg("is_lambda1")=(bool const)(false), bp::arg("map0")=SireBase::PropertyMap(), bp::arg("map1")=SireBase::PropertyMap() )
            , "Update the coordinates and velocities of original_system with those from\nupdated_system.\n\nPar:am system_original\nThe original system.\n\nPar:am system_renumbered\nThe original system, atoms and residues have been renumbered to be\nunique and in ascending order.\n\nPar:am system_updated\nThe updated system, where molecules may not be in the same order.\n\nPar:am map0\nA dictionary of user-defined molecular property names for system0.\n\nPar:am map1\nA dictionary of user-defined molecular property names for system1.\n\nRetval: system, mapping\nThe system with updated coordinates and velocities and a mapping\nbetween the molecule indices in both systems.\n" );
    
    }

}
