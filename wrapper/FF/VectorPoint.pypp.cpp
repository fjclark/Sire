// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "VectorPoint.pypp.hpp"

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

SireFF::VectorPoint __copy__(const SireFF::VectorPoint &other){ return SireFF::VectorPoint(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_VectorPoint_class(){

    { //::SireFF::VectorPoint
        typedef bp::class_< SireFF::VectorPoint, bp::bases< SireFF::Point, SireBase::Property > > VectorPoint_exposer_t;
        VectorPoint_exposer_t VectorPoint_exposer = VectorPoint_exposer_t( "VectorPoint", "This is a simple wrapper for a point in space", bp::init< >("Constructor") );
        bp::scope VectorPoint_scope( VectorPoint_exposer );
        VectorPoint_exposer.def( bp::init< SireMaths::Vector const & >(( bp::arg("point") ), "Constructor for the specified point") );
        VectorPoint_exposer.def( bp::init< SireFF::VectorPoint const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireFF::VectorPoint::addForce
        
            typedef bool ( ::SireFF::VectorPoint::*addForce_function_type)( ::SireFF::MolForceTable &,::SireMaths::Vector const & ) const;
            addForce_function_type addForce_function_value( &::SireFF::VectorPoint::addForce );
            
            VectorPoint_exposer.def( 
                "addForce"
                , addForce_function_value
                , ( bp::arg("molforces"), bp::arg("force") )
                , "No forces on a point" );
        
        }
        { //::SireFF::VectorPoint::addForce
        
            typedef bool ( ::SireFF::VectorPoint::*addForce_function_type)( ::SireFF::ForceTable &,::SireMaths::Vector const & ) const;
            addForce_function_type addForce_function_value( &::SireFF::VectorPoint::addForce );
            
            VectorPoint_exposer.def( 
                "addForce"
                , addForce_function_value
                , ( bp::arg("forces"), bp::arg("force") )
                , "No forces on a point" );
        
        }
        { //::SireFF::VectorPoint::contains
        
            typedef bool ( ::SireFF::VectorPoint::*contains_function_type)( ::SireMol::MolNum ) const;
            contains_function_type contains_function_value( &::SireFF::VectorPoint::contains );
            
            VectorPoint_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molnum") )
                , "No molecules are needed to create this point" );
        
        }
        { //::SireFF::VectorPoint::contains
        
            typedef bool ( ::SireFF::VectorPoint::*contains_function_type)( ::SireMol::MolID const & ) const;
            contains_function_type contains_function_value( &::SireFF::VectorPoint::contains );
            
            VectorPoint_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molid") )
                , "No molecules are needed to create this point" );
        
        }
        { //::SireFF::VectorPoint::isExtraMoleculePoint
        
            typedef bool ( ::SireFF::VectorPoint::*isExtraMoleculePoint_function_type)(  ) const;
            isExtraMoleculePoint_function_type isExtraMoleculePoint_function_value( &::SireFF::VectorPoint::isExtraMoleculePoint );
            
            VectorPoint_exposer.def( 
                "isExtraMoleculePoint"
                , isExtraMoleculePoint_function_value
                , "Return whether or not this is an extramolecular point (it is independent\nof the coordinates of atoms in any molecule, i.e. it is just a point in space)" );
        
        }
        { //::SireFF::VectorPoint::isInterMoleculePoint
        
            typedef bool ( ::SireFF::VectorPoint::*isInterMoleculePoint_function_type)(  ) const;
            isInterMoleculePoint_function_type isInterMoleculePoint_function_value( &::SireFF::VectorPoint::isInterMoleculePoint );
            
            VectorPoint_exposer.def( 
                "isInterMoleculePoint"
                , isInterMoleculePoint_function_value
                , "Return whether or not this is an intermolecular point (it depends on\ncoordinates of atoms from than one molecule)" );
        
        }
        { //::SireFF::VectorPoint::isIntraMoleculePoint
        
            typedef bool ( ::SireFF::VectorPoint::*isIntraMoleculePoint_function_type)(  ) const;
            isIntraMoleculePoint_function_type isIntraMoleculePoint_function_value( &::SireFF::VectorPoint::isIntraMoleculePoint );
            
            VectorPoint_exposer.def( 
                "isIntraMoleculePoint"
                , isIntraMoleculePoint_function_value
                , "Return whether this is an intramolecular point (it depends on coordinates\nof atoms in just one molecule)" );
        
        }
        { //::SireFF::VectorPoint::molecules
        
            typedef ::SireMol::Molecules ( ::SireFF::VectorPoint::*molecules_function_type)(  ) const;
            molecules_function_type molecules_function_value( &::SireFF::VectorPoint::molecules );
            
            VectorPoint_exposer.def( 
                "molecules"
                , molecules_function_value
                , "No molecules are needed to create this point" );
        
        }
        { //::SireFF::VectorPoint::nMolecules
        
            typedef int ( ::SireFF::VectorPoint::*nMolecules_function_type)(  ) const;
            nMolecules_function_type nMolecules_function_value( &::SireFF::VectorPoint::nMolecules );
            
            VectorPoint_exposer.def( 
                "nMolecules"
                , nMolecules_function_value
                , "No molecules are needed to create this point" );
        
        }
        VectorPoint_exposer.def( bp::self != bp::self );
        { //::SireFF::VectorPoint::operator=
        
            typedef ::SireFF::VectorPoint & ( ::SireFF::VectorPoint::*assign_function_type)( ::SireFF::VectorPoint const & ) ;
            assign_function_type assign_function_value( &::SireFF::VectorPoint::operator= );
            
            VectorPoint_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        VectorPoint_exposer.def( bp::self == bp::self );
        { //::SireFF::VectorPoint::toString
        
            typedef ::QString ( ::SireFF::VectorPoint::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireFF::VectorPoint::toString );
            
            VectorPoint_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation" );
        
        }
        { //::SireFF::VectorPoint::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::VectorPoint::typeName );
            
            VectorPoint_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireFF::VectorPoint::update
        
            typedef bool ( ::SireFF::VectorPoint::*update_function_type)( ::SireMol::MoleculeData const & ) ;
            update_function_type update_function_value( &::SireFF::VectorPoint::update );
            
            VectorPoint_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata") )
                , "A VectorPoint is not updatable" );
        
        }
        { //::SireFF::VectorPoint::update
        
            typedef bool ( ::SireFF::VectorPoint::*update_function_type)( ::SireMol::Molecules const & ) ;
            update_function_type update_function_value( &::SireFF::VectorPoint::update );
            
            VectorPoint_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molecules") )
                , "A VectorPoint is not updatable" );
        
        }
        { //::SireFF::VectorPoint::update
        
            typedef bool ( ::SireFF::VectorPoint::*update_function_type)( ::SireMol::MoleculeGroup const & ) ;
            update_function_type update_function_value( &::SireFF::VectorPoint::update );
            
            VectorPoint_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molgroup") )
                , "A VectorPoint is not updatable" );
        
        }
        { //::SireFF::VectorPoint::update
        
            typedef bool ( ::SireFF::VectorPoint::*update_function_type)( ::SireMol::MolGroupsBase const & ) ;
            update_function_type update_function_value( &::SireFF::VectorPoint::update );
            
            VectorPoint_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molgroups") )
                , "A VectorPoint is not updatable" );
        
        }
        { //::SireFF::VectorPoint::usesMoleculesIn
        
            typedef bool ( ::SireFF::VectorPoint::*usesMoleculesIn_function_type)( ::SireFF::ForceTable const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::VectorPoint::usesMoleculesIn );
            
            VectorPoint_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("forcetable") )
                , "No molecules are needed to create this point" );
        
        }
        { //::SireFF::VectorPoint::usesMoleculesIn
        
            typedef bool ( ::SireFF::VectorPoint::*usesMoleculesIn_function_type)( ::SireMol::Molecules const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::VectorPoint::usesMoleculesIn );
            
            VectorPoint_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molecules") )
                , "No molecules are needed to create this point" );
        
        }
        { //::SireFF::VectorPoint::usesMoleculesIn
        
            typedef bool ( ::SireFF::VectorPoint::*usesMoleculesIn_function_type)( ::SireMol::MoleculeGroup const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::VectorPoint::usesMoleculesIn );
            
            VectorPoint_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molgroup") )
                , "No molecules are needed to create this point" );
        
        }
        { //::SireFF::VectorPoint::usesMoleculesIn
        
            typedef bool ( ::SireFF::VectorPoint::*usesMoleculesIn_function_type)( ::SireMol::MolGroupsBase const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::VectorPoint::usesMoleculesIn );
            
            VectorPoint_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molgroups") )
                , "No molecules are needed to create this point" );
        
        }
        { //::SireFF::VectorPoint::wouldUpdate
        
            typedef bool ( ::SireFF::VectorPoint::*wouldUpdate_function_type)( ::SireMol::MoleculeData const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::VectorPoint::wouldUpdate );
            
            VectorPoint_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("moldata") )
                , "A VectorPoint is not updatable" );
        
        }
        { //::SireFF::VectorPoint::wouldUpdate
        
            typedef bool ( ::SireFF::VectorPoint::*wouldUpdate_function_type)( ::SireMol::Molecules const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::VectorPoint::wouldUpdate );
            
            VectorPoint_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molecules") )
                , "A VectorPoint is not updatable" );
        
        }
        { //::SireFF::VectorPoint::wouldUpdate
        
            typedef bool ( ::SireFF::VectorPoint::*wouldUpdate_function_type)( ::SireMol::MoleculeGroup const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::VectorPoint::wouldUpdate );
            
            VectorPoint_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molgroup") )
                , "A VectorPoint is not updatable" );
        
        }
        { //::SireFF::VectorPoint::wouldUpdate
        
            typedef bool ( ::SireFF::VectorPoint::*wouldUpdate_function_type)( ::SireMol::MolGroupsBase const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::VectorPoint::wouldUpdate );
            
            VectorPoint_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molgroups") )
                , "A VectorPoint is not updatable" );
        
        }
        VectorPoint_exposer.staticmethod( "typeName" );
        VectorPoint_exposer.def( "__copy__", &__copy__);
        VectorPoint_exposer.def( "__deepcopy__", &__copy__);
        VectorPoint_exposer.def( "clone", &__copy__);
        VectorPoint_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::VectorPoint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        VectorPoint_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::VectorPoint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        VectorPoint_exposer.def( "__setstate__", &__setstate__base64< ::SireFF::VectorPoint > );
        VectorPoint_exposer.def( "__getstate__", &__getstate__base64< ::SireFF::VectorPoint > );
        VectorPoint_exposer.def( "__str__", &__str__< ::SireFF::VectorPoint > );
        VectorPoint_exposer.def( "__repr__", &__str__< ::SireFF::VectorPoint > );
    }

}
