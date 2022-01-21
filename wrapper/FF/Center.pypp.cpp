// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Center.pypp.hpp"

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

SireFF::Center __copy__(const SireFF::Center &other){ return SireFF::Center(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Center_class(){

    { //::SireFF::Center
        typedef bp::class_< SireFF::Center, bp::bases< SireFF::Point, SireBase::Property > > Center_exposer_t;
        Center_exposer_t Center_exposer = Center_exposer_t( "Center", "This point returns the center of a view of a molecule, or group\nof molecules", bp::init< >("Constructor") );
        bp::scope Center_scope( Center_exposer );
        Center_exposer.def( bp::init< SireMol::MoleculeView const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molview"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to get the center of the molecule view molview using the\npassed property map to find the required properties\nThrow: SireBase::missing_property\nThrow: SireError::invalid_cast\n") );
        Center_exposer.def( bp::init< SireMol::Molecules const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecules"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to get the center of the molecules in molecules, using the\npassed property map to find the required properties\nThrow: SireBase::missing_property\nThrow: SireError::invalid_cast\n") );
        Center_exposer.def( bp::init< SireFF::Center const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireFF::Center::addForce
        
            typedef bool ( ::SireFF::Center::*addForce_function_type)( ::SireFF::MolForceTable &,::SireMaths::Vector const & ) const;
            addForce_function_type addForce_function_value( &::SireFF::Center::addForce );
            
            Center_exposer.def( 
                "addForce"
                , addForce_function_value
                , ( bp::arg("molforces"), bp::arg("force") )
                , "Decompose the force force acting on this point from the\nmolecule whose forces are in molforces and add the\nforce onto the table" );
        
        }
        { //::SireFF::Center::addForce
        
            typedef bool ( ::SireFF::Center::*addForce_function_type)( ::SireFF::ForceTable &,::SireMaths::Vector const & ) const;
            addForce_function_type addForce_function_value( &::SireFF::Center::addForce );
            
            Center_exposer.def( 
                "addForce"
                , addForce_function_value
                , ( bp::arg("forces"), bp::arg("force") )
                , "Decompose the force force into the forces acting on\nthe molecules that contribute to this point and add those\nforces onto the table forces" );
        
        }
        { //::SireFF::Center::contains
        
            typedef bool ( ::SireFF::Center::*contains_function_type)( ::SireMol::MolNum ) const;
            contains_function_type contains_function_value( &::SireFF::Center::contains );
            
            Center_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molnum") )
                , "Return whether or not the molecule with number molnum is\nneeded to generate this point" );
        
        }
        { //::SireFF::Center::contains
        
            typedef bool ( ::SireFF::Center::*contains_function_type)( ::SireMol::MolID const & ) const;
            contains_function_type contains_function_value( &::SireFF::Center::contains );
            
            Center_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molid") )
                , "Return whether or not this molecule with ID molid is\nneeded to generate this point" );
        
        }
        { //::SireFF::Center::isExtraMoleculePoint
        
            typedef bool ( ::SireFF::Center::*isExtraMoleculePoint_function_type)(  ) const;
            isExtraMoleculePoint_function_type isExtraMoleculePoint_function_value( &::SireFF::Center::isExtraMoleculePoint );
            
            Center_exposer.def( 
                "isExtraMoleculePoint"
                , isExtraMoleculePoint_function_value
                , "Return whether or not this is an extramolecular point (it is independent\nof the coordinates of atoms in any molecule, i.e. it is just a point in space)" );
        
        }
        { //::SireFF::Center::isInterMoleculePoint
        
            typedef bool ( ::SireFF::Center::*isInterMoleculePoint_function_type)(  ) const;
            isInterMoleculePoint_function_type isInterMoleculePoint_function_value( &::SireFF::Center::isInterMoleculePoint );
            
            Center_exposer.def( 
                "isInterMoleculePoint"
                , isInterMoleculePoint_function_value
                , "Return whether or not this is an intermolecular point (it depends on\ncoordinates of atoms from than one molecule)" );
        
        }
        { //::SireFF::Center::isIntraMoleculePoint
        
            typedef bool ( ::SireFF::Center::*isIntraMoleculePoint_function_type)(  ) const;
            isIntraMoleculePoint_function_type isIntraMoleculePoint_function_value( &::SireFF::Center::isIntraMoleculePoint );
            
            Center_exposer.def( 
                "isIntraMoleculePoint"
                , isIntraMoleculePoint_function_value
                , "Return whether this is an intramolecular point (it depends on coordinates\nof atoms in just one molecule)" );
        
        }
        { //::SireFF::Center::molecules
        
            typedef ::SireMol::Molecules ( ::SireFF::Center::*molecules_function_type)(  ) const;
            molecules_function_type molecules_function_value( &::SireFF::Center::molecules );
            
            Center_exposer.def( 
                "molecules"
                , molecules_function_value
                , "Return all of the molecules used to generate this point" );
        
        }
        { //::SireFF::Center::nMolecules
        
            typedef int ( ::SireFF::Center::*nMolecules_function_type)(  ) const;
            nMolecules_function_type nMolecules_function_value( &::SireFF::Center::nMolecules );
            
            Center_exposer.def( 
                "nMolecules"
                , nMolecules_function_value
                , "Return the number of molecules needed to generate this point" );
        
        }
        Center_exposer.def( bp::self != bp::self );
        { //::SireFF::Center::operator=
        
            typedef ::SireFF::Center & ( ::SireFF::Center::*assign_function_type)( ::SireFF::Center const & ) ;
            assign_function_type assign_function_value( &::SireFF::Center::operator= );
            
            Center_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Center_exposer.def( bp::self == bp::self );
        { //::SireFF::Center::setSpace
        
            typedef void ( ::SireFF::Center::*setSpace_function_type)( ::SireVol::Space const & ) ;
            setSpace_function_type setSpace_function_value( &::SireFF::Center::setSpace );
            
            Center_exposer.def( 
                "setSpace"
                , setSpace_function_value
                , ( bp::arg("space") )
                , "Set the space - if there is more than one molecule, then this\npoint can only be used with a cartesian, non-periodic space" );
        
        }
        { //::SireFF::Center::toString
        
            typedef ::QString ( ::SireFF::Center::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireFF::Center::toString );
            
            Center_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation" );
        
        }
        { //::SireFF::Center::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Center::typeName );
            
            Center_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireFF::Center::update
        
            typedef bool ( ::SireFF::Center::*update_function_type)( ::SireMol::MoleculeData const & ) ;
            update_function_type update_function_value( &::SireFF::Center::update );
            
            Center_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata") )
                , "Update the molecules used to create this point" );
        
        }
        { //::SireFF::Center::update
        
            typedef bool ( ::SireFF::Center::*update_function_type)( ::SireMol::Molecules const & ) ;
            update_function_type update_function_value( &::SireFF::Center::update );
            
            Center_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molecules") )
                , "Update the molecules used to create this point" );
        
        }
        { //::SireFF::Center::update
        
            typedef bool ( ::SireFF::Center::*update_function_type)( ::SireMol::MoleculeGroup const & ) ;
            update_function_type update_function_value( &::SireFF::Center::update );
            
            Center_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molgroup") )
                , "Update the molecules used to create this point" );
        
        }
        { //::SireFF::Center::update
        
            typedef bool ( ::SireFF::Center::*update_function_type)( ::SireMol::MolGroupsBase const & ) ;
            update_function_type update_function_value( &::SireFF::Center::update );
            
            Center_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molgroups") )
                , "Update the molecules used to create this point" );
        
        }
        { //::SireFF::Center::usesMoleculesIn
        
            typedef bool ( ::SireFF::Center::*usesMoleculesIn_function_type)( ::SireFF::ForceTable const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::Center::usesMoleculesIn );
            
            Center_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("forcetable") )
                , "Return whether or not this point uses data from any of the\nmolecules in the passed forcetable" );
        
        }
        { //::SireFF::Center::usesMoleculesIn
        
            typedef bool ( ::SireFF::Center::*usesMoleculesIn_function_type)( ::SireMol::Molecules const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::Center::usesMoleculesIn );
            
            Center_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molecules") )
                , "Return whether or not this point uses data from any of the\nmolecules in molecules" );
        
        }
        { //::SireFF::Center::usesMoleculesIn
        
            typedef bool ( ::SireFF::Center::*usesMoleculesIn_function_type)( ::SireMol::MoleculeGroup const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::Center::usesMoleculesIn );
            
            Center_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molgroup") )
                , "Return whether or not this point uses data from any of the\nmolecules in the group molgroup" );
        
        }
        { //::SireFF::Center::usesMoleculesIn
        
            typedef bool ( ::SireFF::Center::*usesMoleculesIn_function_type)( ::SireMol::MolGroupsBase const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::Center::usesMoleculesIn );
            
            Center_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molgroups") )
                , "Return whether or not this point uses data from any of the\nmolecules in the groups in molgroups" );
        
        }
        { //::SireFF::Center::wouldUpdate
        
            typedef bool ( ::SireFF::Center::*wouldUpdate_function_type)( ::SireMol::MoleculeData const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::Center::wouldUpdate );
            
            Center_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("moldata") )
                , "Return whether or not the passed molecule would change this point" );
        
        }
        { //::SireFF::Center::wouldUpdate
        
            typedef bool ( ::SireFF::Center::*wouldUpdate_function_type)( ::SireMol::Molecules const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::Center::wouldUpdate );
            
            Center_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molecules") )
                , "Return whether or not the passed molecules would change this point" );
        
        }
        { //::SireFF::Center::wouldUpdate
        
            typedef bool ( ::SireFF::Center::*wouldUpdate_function_type)( ::SireMol::MoleculeGroup const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::Center::wouldUpdate );
            
            Center_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molgroup") )
                , "Return whether or not the passed molecules would change this point" );
        
        }
        { //::SireFF::Center::wouldUpdate
        
            typedef bool ( ::SireFF::Center::*wouldUpdate_function_type)( ::SireMol::MolGroupsBase const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::Center::wouldUpdate );
            
            Center_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molgroups") )
                , "Return whether or not the passed molecules would change this point" );
        
        }
        Center_exposer.staticmethod( "typeName" );
        Center_exposer.def( "__copy__", &__copy__);
        Center_exposer.def( "__deepcopy__", &__copy__);
        Center_exposer.def( "clone", &__copy__);
        Center_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::Center >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Center_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::Center >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Center_exposer.def( "__getstate_manages_dict__", true);
        Center_exposer.def( "__safe_for_unpickling__", true);
        Center_exposer.def( "__setstate__", &__setstate__base64< ::SireFF::Center > );
        Center_exposer.def( "__getstate__", &__getstate__base64< ::SireFF::Center > );
        Center_exposer.def( "__str__", &__str__< ::SireFF::Center > );
        Center_exposer.def( "__repr__", &__str__< ::SireFF::Center > );
    }

}
