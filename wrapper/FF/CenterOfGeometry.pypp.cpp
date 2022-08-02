// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CenterOfGeometry.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/evaluator.h"

#include "SireMol/mgidx.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/moleculegroups.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireVol/aabox.h"

#include "SireVol/errors.h"

#include "forcetable.h"

#include "point.h"

#include <boost/tuple/tuple.hpp>

#include "point.h"

SireFF::CenterOfGeometry __copy__(const SireFF::CenterOfGeometry &other){ return SireFF::CenterOfGeometry(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_CenterOfGeometry_class(){

    { //::SireFF::CenterOfGeometry
        typedef bp::class_< SireFF::CenterOfGeometry, bp::bases< SireFF::Point, SireBase::Property > > CenterOfGeometry_exposer_t;
        CenterOfGeometry_exposer_t CenterOfGeometry_exposer = CenterOfGeometry_exposer_t( "CenterOfGeometry", "This point returns the center of geometry of a view of a molecule,\nor group of molecules", bp::init< >("Constructor") );
        bp::scope CenterOfGeometry_scope( CenterOfGeometry_exposer );
        CenterOfGeometry_exposer.def( bp::init< SireMol::MoleculeView const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molview"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to get the center of the molecule view molview using the\npassed property map to find the required properties\nThrow: SireBase::missing_property\nThrow: SireError::invalid_cast\n") );
        CenterOfGeometry_exposer.def( bp::init< SireMol::Molecules const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecules"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to get the center of the molecules in molecules, using the\npassed property map to find the required properties\nThrow: SireBase::missing_property\nThrow: SireError::invalid_cast\n") );
        CenterOfGeometry_exposer.def( bp::init< SireFF::CenterOfGeometry const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireFF::CenterOfGeometry::addForce
        
            typedef bool ( ::SireFF::CenterOfGeometry::*addForce_function_type)( ::SireFF::MolForceTable &,::SireMaths::Vector const & ) const;
            addForce_function_type addForce_function_value( &::SireFF::CenterOfGeometry::addForce );
            
            CenterOfGeometry_exposer.def( 
                "addForce"
                , addForce_function_value
                , ( bp::arg("molforces"), bp::arg("force") )
                , bp::release_gil_policy()
                , "Decompose the force force acting on this point from the\nmolecule whose forces are in molforces and add the\nforce onto the table" );
        
        }
        { //::SireFF::CenterOfGeometry::addForce
        
            typedef bool ( ::SireFF::CenterOfGeometry::*addForce_function_type)( ::SireFF::ForceTable &,::SireMaths::Vector const & ) const;
            addForce_function_type addForce_function_value( &::SireFF::CenterOfGeometry::addForce );
            
            CenterOfGeometry_exposer.def( 
                "addForce"
                , addForce_function_value
                , ( bp::arg("forces"), bp::arg("force") )
                , bp::release_gil_policy()
                , "Decompose the force force into the forces acting on\nthe molecules that contribute to this point and add those\nforces onto the table forces" );
        
        }
        { //::SireFF::CenterOfGeometry::contains
        
            typedef bool ( ::SireFF::CenterOfGeometry::*contains_function_type)( ::SireMol::MolNum ) const;
            contains_function_type contains_function_value( &::SireFF::CenterOfGeometry::contains );
            
            CenterOfGeometry_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molnum") )
                , bp::release_gil_policy()
                , "Return whether or not the molecule with number molnum is\nneeded to generate this point" );
        
        }
        { //::SireFF::CenterOfGeometry::contains
        
            typedef bool ( ::SireFF::CenterOfGeometry::*contains_function_type)( ::SireMol::MolID const & ) const;
            contains_function_type contains_function_value( &::SireFF::CenterOfGeometry::contains );
            
            CenterOfGeometry_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molid") )
                , bp::release_gil_policy()
                , "Return whether or not this molecule with ID molid is\nneeded to generate this point" );
        
        }
        { //::SireFF::CenterOfGeometry::isExtraMoleculePoint
        
            typedef bool ( ::SireFF::CenterOfGeometry::*isExtraMoleculePoint_function_type)(  ) const;
            isExtraMoleculePoint_function_type isExtraMoleculePoint_function_value( &::SireFF::CenterOfGeometry::isExtraMoleculePoint );
            
            CenterOfGeometry_exposer.def( 
                "isExtraMoleculePoint"
                , isExtraMoleculePoint_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is an extramolecular point (it is independent\nof the coordinates of atoms in any molecule, i.e. it is just a point in space)" );
        
        }
        { //::SireFF::CenterOfGeometry::isInterMoleculePoint
        
            typedef bool ( ::SireFF::CenterOfGeometry::*isInterMoleculePoint_function_type)(  ) const;
            isInterMoleculePoint_function_type isInterMoleculePoint_function_value( &::SireFF::CenterOfGeometry::isInterMoleculePoint );
            
            CenterOfGeometry_exposer.def( 
                "isInterMoleculePoint"
                , isInterMoleculePoint_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is an intermolecular point (it depends on\ncoordinates of atoms from than one molecule)" );
        
        }
        { //::SireFF::CenterOfGeometry::isIntraMoleculePoint
        
            typedef bool ( ::SireFF::CenterOfGeometry::*isIntraMoleculePoint_function_type)(  ) const;
            isIntraMoleculePoint_function_type isIntraMoleculePoint_function_value( &::SireFF::CenterOfGeometry::isIntraMoleculePoint );
            
            CenterOfGeometry_exposer.def( 
                "isIntraMoleculePoint"
                , isIntraMoleculePoint_function_value
                , bp::release_gil_policy()
                , "Return whether this is an intramolecular point (it depends on coordinates\nof atoms in just one molecule)" );
        
        }
        { //::SireFF::CenterOfGeometry::molecules
        
            typedef ::SireMol::Molecules ( ::SireFF::CenterOfGeometry::*molecules_function_type)(  ) const;
            molecules_function_type molecules_function_value( &::SireFF::CenterOfGeometry::molecules );
            
            CenterOfGeometry_exposer.def( 
                "molecules"
                , molecules_function_value
                , bp::release_gil_policy()
                , "Return all of the molecules used to generate this point" );
        
        }
        { //::SireFF::CenterOfGeometry::nMolecules
        
            typedef int ( ::SireFF::CenterOfGeometry::*nMolecules_function_type)(  ) const;
            nMolecules_function_type nMolecules_function_value( &::SireFF::CenterOfGeometry::nMolecules );
            
            CenterOfGeometry_exposer.def( 
                "nMolecules"
                , nMolecules_function_value
                , bp::release_gil_policy()
                , "Return the number of molecules needed to generate this point" );
        
        }
        CenterOfGeometry_exposer.def( bp::self != bp::self );
        { //::SireFF::CenterOfGeometry::operator=
        
            typedef ::SireFF::CenterOfGeometry & ( ::SireFF::CenterOfGeometry::*assign_function_type)( ::SireFF::CenterOfGeometry const & ) ;
            assign_function_type assign_function_value( &::SireFF::CenterOfGeometry::operator= );
            
            CenterOfGeometry_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CenterOfGeometry_exposer.def( bp::self == bp::self );
        { //::SireFF::CenterOfGeometry::setSpace
        
            typedef void ( ::SireFF::CenterOfGeometry::*setSpace_function_type)( ::SireVol::Space const & ) ;
            setSpace_function_type setSpace_function_value( &::SireFF::CenterOfGeometry::setSpace );
            
            CenterOfGeometry_exposer.def( 
                "setSpace"
                , setSpace_function_value
                , ( bp::arg("space") )
                , bp::release_gil_policy()
                , "Set the space used by this point - a CenterOfGeometry cannot\nbe calculated for periodic or non-cartesian spaces if there\nis more than one molecule" );
        
        }
        { //::SireFF::CenterOfGeometry::toString
        
            typedef ::QString ( ::SireFF::CenterOfGeometry::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireFF::CenterOfGeometry::toString );
            
            CenterOfGeometry_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation" );
        
        }
        { //::SireFF::CenterOfGeometry::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::CenterOfGeometry::typeName );
            
            CenterOfGeometry_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::CenterOfGeometry::update
        
            typedef bool ( ::SireFF::CenterOfGeometry::*update_function_type)( ::SireMol::MoleculeData const & ) ;
            update_function_type update_function_value( &::SireFF::CenterOfGeometry::update );
            
            CenterOfGeometry_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata") )
                , bp::release_gil_policy()
                , "Update the molecules used to create this point" );
        
        }
        { //::SireFF::CenterOfGeometry::update
        
            typedef bool ( ::SireFF::CenterOfGeometry::*update_function_type)( ::SireMol::Molecules const & ) ;
            update_function_type update_function_value( &::SireFF::CenterOfGeometry::update );
            
            CenterOfGeometry_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molecules") )
                , bp::release_gil_policy()
                , "Update the molecules used to create this point" );
        
        }
        { //::SireFF::CenterOfGeometry::update
        
            typedef bool ( ::SireFF::CenterOfGeometry::*update_function_type)( ::SireMol::MoleculeGroup const & ) ;
            update_function_type update_function_value( &::SireFF::CenterOfGeometry::update );
            
            CenterOfGeometry_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molgroup") )
                , bp::release_gil_policy()
                , "Update the molecules used to create this point" );
        
        }
        { //::SireFF::CenterOfGeometry::update
        
            typedef bool ( ::SireFF::CenterOfGeometry::*update_function_type)( ::SireMol::MolGroupsBase const & ) ;
            update_function_type update_function_value( &::SireFF::CenterOfGeometry::update );
            
            CenterOfGeometry_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molgroups") )
                , bp::release_gil_policy()
                , "Update the molecules used to create this point" );
        
        }
        { //::SireFF::CenterOfGeometry::usesMoleculesIn
        
            typedef bool ( ::SireFF::CenterOfGeometry::*usesMoleculesIn_function_type)( ::SireFF::ForceTable const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::CenterOfGeometry::usesMoleculesIn );
            
            CenterOfGeometry_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("forcetable") )
                , bp::release_gil_policy()
                , "Return whether or not this point uses data from any of the\nmolecules in the passed forcetable" );
        
        }
        { //::SireFF::CenterOfGeometry::usesMoleculesIn
        
            typedef bool ( ::SireFF::CenterOfGeometry::*usesMoleculesIn_function_type)( ::SireMol::Molecules const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::CenterOfGeometry::usesMoleculesIn );
            
            CenterOfGeometry_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molecules") )
                , bp::release_gil_policy()
                , "Return whether or not this point uses data from any of the\nmolecules in molecules" );
        
        }
        { //::SireFF::CenterOfGeometry::usesMoleculesIn
        
            typedef bool ( ::SireFF::CenterOfGeometry::*usesMoleculesIn_function_type)( ::SireMol::MoleculeGroup const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::CenterOfGeometry::usesMoleculesIn );
            
            CenterOfGeometry_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molgroup") )
                , bp::release_gil_policy()
                , "Return whether or not this point uses data from any of the\nmolecules in the group molgroup" );
        
        }
        { //::SireFF::CenterOfGeometry::usesMoleculesIn
        
            typedef bool ( ::SireFF::CenterOfGeometry::*usesMoleculesIn_function_type)( ::SireMol::MolGroupsBase const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::CenterOfGeometry::usesMoleculesIn );
            
            CenterOfGeometry_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molgroups") )
                , bp::release_gil_policy()
                , "Return whether or not this point uses data from any of the\nmolecules in the groups in molgroups" );
        
        }
        { //::SireFF::CenterOfGeometry::wouldUpdate
        
            typedef bool ( ::SireFF::CenterOfGeometry::*wouldUpdate_function_type)( ::SireMol::MoleculeData const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::CenterOfGeometry::wouldUpdate );
            
            CenterOfGeometry_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("moldata") )
                , bp::release_gil_policy()
                , "Return whether or not the passed molecule would change this point" );
        
        }
        { //::SireFF::CenterOfGeometry::wouldUpdate
        
            typedef bool ( ::SireFF::CenterOfGeometry::*wouldUpdate_function_type)( ::SireMol::Molecules const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::CenterOfGeometry::wouldUpdate );
            
            CenterOfGeometry_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molecules") )
                , bp::release_gil_policy()
                , "Return whether or not the passed molecules would change this point" );
        
        }
        { //::SireFF::CenterOfGeometry::wouldUpdate
        
            typedef bool ( ::SireFF::CenterOfGeometry::*wouldUpdate_function_type)( ::SireMol::MoleculeGroup const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::CenterOfGeometry::wouldUpdate );
            
            CenterOfGeometry_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molgroup") )
                , bp::release_gil_policy()
                , "Return whether or not the passed molecules would change this point" );
        
        }
        { //::SireFF::CenterOfGeometry::wouldUpdate
        
            typedef bool ( ::SireFF::CenterOfGeometry::*wouldUpdate_function_type)( ::SireMol::MolGroupsBase const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::CenterOfGeometry::wouldUpdate );
            
            CenterOfGeometry_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molgroups") )
                , bp::release_gil_policy()
                , "Return whether or not the passed molecules would change this point" );
        
        }
        CenterOfGeometry_exposer.staticmethod( "typeName" );
        CenterOfGeometry_exposer.def( "__copy__", &__copy__);
        CenterOfGeometry_exposer.def( "__deepcopy__", &__copy__);
        CenterOfGeometry_exposer.def( "clone", &__copy__);
        CenterOfGeometry_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::CenterOfGeometry >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CenterOfGeometry_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::CenterOfGeometry >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CenterOfGeometry_exposer.def_pickle(sire_pickle_suite< ::SireFF::CenterOfGeometry >());
        CenterOfGeometry_exposer.def( "__str__", &__str__< ::SireFF::CenterOfGeometry > );
        CenterOfGeometry_exposer.def( "__repr__", &__str__< ::SireFF::CenterOfGeometry > );
    }

}
