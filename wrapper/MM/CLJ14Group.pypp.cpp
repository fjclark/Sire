// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "CLJ14Group.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMM/atomljs.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/core.h"

#include "SireMol/editor.hpp"

#include "SireMol/mover.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "clj14group.h"

#include "cljnbpairs.h"

#include "clj14group.h"

SireMM::CLJ14Group __copy__(const SireMM::CLJ14Group &other){ return SireMM::CLJ14Group(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_CLJ14Group_class(){

    { //::SireMM::CLJ14Group
        typedef bp::class_< SireMM::CLJ14Group > CLJ14Group_exposer_t;
        CLJ14Group_exposer_t CLJ14Group_exposer = CLJ14Group_exposer_t( "CLJ14Group", "This class holds all of the information needed to calculate\nthe 14-nonbonded energy for a molecule\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope CLJ14Group_scope( CLJ14Group_exposer );
        CLJ14Group_exposer.def( bp::init< SireMol::MoleculeView const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to calculate the 14-energy of the passed molecule") );
        CLJ14Group_exposer.def( bp::init< SireMol::MoleculeView const &, SireMM::CLJFunction::COMBINING_RULES, bool, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("combining_rules"), bp::arg("is_strict"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to calculate the 14-energy of the passed molecule using the\nsupplied combining rules and strict mode") );
        CLJ14Group_exposer.def( bp::init< SireMM::CLJ14Group const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::CLJ14Group::add
        
            typedef void ( ::SireMM::CLJ14Group::*add_function_type)( ::SireMol::MoleculeView const & ) ;
            add_function_type add_function_value( &::SireMM::CLJ14Group::add );
            
            CLJ14Group_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("new_molecule") )
                , bp::release_gil_policy()
                , "Add the passed molecule to this group" );
        
        }
        { //::SireMM::CLJ14Group::add
        
            typedef void ( ::SireMM::CLJ14Group::*add_function_type)( ::SireMol::AtomSelection const & ) ;
            add_function_type add_function_value( &::SireMM::CLJ14Group::add );
            
            CLJ14Group_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("new_selection") )
                , bp::release_gil_policy()
                , "Add the passed selection onto this group" );
        
        }
        { //::SireMM::CLJ14Group::combiningRules
        
            typedef ::SireMM::CLJFunction::COMBINING_RULES ( ::SireMM::CLJ14Group::*combiningRules_function_type)(  ) const;
            combiningRules_function_type combiningRules_function_value( &::SireMM::CLJ14Group::combiningRules );
            
            CLJ14Group_exposer.def( 
                "combiningRules"
                , combiningRules_function_value
                , bp::release_gil_policy()
                , "Return the type of combining rules in place" );
        
        }
        { //::SireMM::CLJ14Group::energy
        
            typedef ::boost::tuples::tuple< double, double, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( ::SireMM::CLJ14Group::*energy_function_type)(  ) ;
            energy_function_type energy_function_value( &::SireMM::CLJ14Group::energy );
            
            CLJ14Group_exposer.def( 
                "energy"
                , energy_function_value
                , bp::release_gil_policy()
                , "Calculate and return the coulomb and LJ 14 energy" );
        
        }
        { //::SireMM::CLJ14Group::isNull
        
            typedef bool ( ::SireMM::CLJ14Group::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMM::CLJ14Group::isNull );
            
            CLJ14Group_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJ14Group::isStrict
        
            typedef bool ( ::SireMM::CLJ14Group::*isStrict_function_type)(  ) const;
            isStrict_function_type isStrict_function_value( &::SireMM::CLJ14Group::isStrict );
            
            CLJ14Group_exposer.def( 
                "isStrict"
                , isStrict_function_value
                , bp::release_gil_policy()
                , "Return whether or not strict mode is on" );
        
        }
        { //::SireMM::CLJ14Group::molecule
        
            typedef ::SireMol::MoleculeView const & ( ::SireMM::CLJ14Group::*molecule_function_type)(  ) const;
            molecule_function_type molecule_function_value( &::SireMM::CLJ14Group::molecule );
            
            CLJ14Group_exposer.def( 
                "molecule"
                , molecule_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the molecule that is in this group" );
        
        }
        { //::SireMM::CLJ14Group::mustNowRecalculateFromScratch
        
            typedef void ( ::SireMM::CLJ14Group::*mustNowRecalculateFromScratch_function_type)(  ) ;
            mustNowRecalculateFromScratch_function_type mustNowRecalculateFromScratch_function_value( &::SireMM::CLJ14Group::mustNowRecalculateFromScratch );
            
            CLJ14Group_exposer.def( 
                "mustNowRecalculateFromScratch"
                , mustNowRecalculateFromScratch_function_value
                , bp::release_gil_policy()
                , "Set the flag to say that we must recalculate everything\nfrom scratch" );
        
        }
        { //::SireMM::CLJ14Group::mustReallyRecalculateFromScratch
        
            typedef void ( ::SireMM::CLJ14Group::*mustReallyRecalculateFromScratch_function_type)(  ) ;
            mustReallyRecalculateFromScratch_function_type mustReallyRecalculateFromScratch_function_value( &::SireMM::CLJ14Group::mustReallyRecalculateFromScratch );
            
            CLJ14Group_exposer.def( 
                "mustReallyRecalculateFromScratch"
                , mustReallyRecalculateFromScratch_function_value
                , bp::release_gil_policy()
                , "Set the flag to ensure that the energy is really completely\nrecalculated from scratch" );
        
        }
        CLJ14Group_exposer.def( bp::self != bp::self );
        { //::SireMM::CLJ14Group::operator=
        
            typedef ::SireMM::CLJ14Group & ( ::SireMM::CLJ14Group::*assign_function_type)( ::SireMM::CLJ14Group const & ) ;
            assign_function_type assign_function_value( &::SireMM::CLJ14Group::operator= );
            
            CLJ14Group_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CLJ14Group_exposer.def( bp::self == bp::self );
        { //::SireMM::CLJ14Group::propertyMap
        
            typedef ::SireBase::PropertyMap ( ::SireMM::CLJ14Group::*propertyMap_function_type)(  ) const;
            propertyMap_function_type propertyMap_function_value( &::SireMM::CLJ14Group::propertyMap );
            
            CLJ14Group_exposer.def( 
                "propertyMap"
                , propertyMap_function_value
                , bp::release_gil_policy()
                , "Return the property map used to find the properties needed to\ncalculate the 14 energy (coordinates, charge, LJ\nand intrascale)" );
        
        }
        { //::SireMM::CLJ14Group::recalculatingFromScratch
        
            typedef bool ( ::SireMM::CLJ14Group::*recalculatingFromScratch_function_type)(  ) const;
            recalculatingFromScratch_function_type recalculatingFromScratch_function_value( &::SireMM::CLJ14Group::recalculatingFromScratch );
            
            CLJ14Group_exposer.def( 
                "recalculatingFromScratch"
                , recalculatingFromScratch_function_value
                , bp::release_gil_policy()
                , "Return whether or not we are recalculating things from scratch" );
        
        }
        { //::SireMM::CLJ14Group::remove
        
            typedef void ( ::SireMM::CLJ14Group::*remove_function_type)( ::SireMol::AtomSelection const & ) ;
            remove_function_type remove_function_value( &::SireMM::CLJ14Group::remove );
            
            CLJ14Group_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("new_selection") )
                , bp::release_gil_policy()
                , "Remove the passed selection from the group" );
        
        }
        { //::SireMM::CLJ14Group::remove
        
            typedef void ( ::SireMM::CLJ14Group::*remove_function_type)( ::SireMol::MoleculeView const & ) ;
            remove_function_type remove_function_value( &::SireMM::CLJ14Group::remove );
            
            CLJ14Group_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("new_molecule") )
                , bp::release_gil_policy()
                , "Remove the passed molecule from this group" );
        
        }
        { //::SireMM::CLJ14Group::setArithmeticCombiningRules
        
            typedef void ( ::SireMM::CLJ14Group::*setArithmeticCombiningRules_function_type)( bool ) ;
            setArithmeticCombiningRules_function_type setArithmeticCombiningRules_function_value( &::SireMM::CLJ14Group::setArithmeticCombiningRules );
            
            CLJ14Group_exposer.def( 
                "setArithmeticCombiningRules"
                , setArithmeticCombiningRules_function_value
                , ( bp::arg("on") )
                , bp::release_gil_policy()
                , "Switch on or off the use of arithmetic combining rules" );
        
        }
        { //::SireMM::CLJ14Group::setCombiningRules
        
            typedef void ( ::SireMM::CLJ14Group::*setCombiningRules_function_type)( ::SireMM::CLJFunction::COMBINING_RULES ) ;
            setCombiningRules_function_type setCombiningRules_function_value( &::SireMM::CLJ14Group::setCombiningRules );
            
            CLJ14Group_exposer.def( 
                "setCombiningRules"
                , setCombiningRules_function_value
                , ( bp::arg("rules") )
                , bp::release_gil_policy()
                , "Set the combining rules to rules" );
        
        }
        { //::SireMM::CLJ14Group::setGeometricCombiningRules
        
            typedef void ( ::SireMM::CLJ14Group::*setGeometricCombiningRules_function_type)( bool ) ;
            setGeometricCombiningRules_function_type setGeometricCombiningRules_function_value( &::SireMM::CLJ14Group::setGeometricCombiningRules );
            
            CLJ14Group_exposer.def( 
                "setGeometricCombiningRules"
                , setGeometricCombiningRules_function_value
                , ( bp::arg("on") )
                , bp::release_gil_policy()
                , "Switch on or off the use og geometric combining rules" );
        
        }
        { //::SireMM::CLJ14Group::setStrict
        
            typedef bool ( ::SireMM::CLJ14Group::*setStrict_function_type)( bool ) ;
            setStrict_function_type setStrict_function_value( &::SireMM::CLJ14Group::setStrict );
            
            CLJ14Group_exposer.def( 
                "setStrict"
                , setStrict_function_value
                , ( bp::arg("isstrict") )
                , bp::release_gil_policy()
                , "Set whether or not strict mode is on. If strict mode is on,\nthen this means that the 1-4 energy is calculated only if both of the\natoms are selected. If strict mode is off, then the 1-4 energy\nis calculated when at least one of the atoms is selected." );
        
        }
        { //::SireMM::CLJ14Group::toString
        
            typedef ::QString ( ::SireMM::CLJ14Group::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::CLJ14Group::toString );
            
            CLJ14Group_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJ14Group::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::CLJ14Group::typeName );
            
            CLJ14Group_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJ14Group::update
        
            typedef void ( ::SireMM::CLJ14Group::*update_function_type)( ::SireMol::MoleculeView const & ) ;
            update_function_type update_function_value( &::SireMM::CLJ14Group::update );
            
            CLJ14Group_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("new_molecule") )
                , bp::release_gil_policy()
                , "Update the contained molecule to the newest version" );
        
        }
        { //::SireMM::CLJ14Group::updateSelection
        
            typedef void ( ::SireMM::CLJ14Group::*updateSelection_function_type)( ::SireMol::AtomSelection const & ) ;
            updateSelection_function_type updateSelection_function_value( &::SireMM::CLJ14Group::updateSelection );
            
            CLJ14Group_exposer.def( 
                "updateSelection"
                , updateSelection_function_value
                , ( bp::arg("selection") )
                , bp::release_gil_policy()
                , "Update the selection to the new passed value" );
        
        }
        { //::SireMM::CLJ14Group::usingArithmeticCombiningRules
        
            typedef bool ( ::SireMM::CLJ14Group::*usingArithmeticCombiningRules_function_type)(  ) const;
            usingArithmeticCombiningRules_function_type usingArithmeticCombiningRules_function_value( &::SireMM::CLJ14Group::usingArithmeticCombiningRules );
            
            CLJ14Group_exposer.def( 
                "usingArithmeticCombiningRules"
                , usingArithmeticCombiningRules_function_value
                , bp::release_gil_policy()
                , "Return whether or not arithmetic combining rules are used" );
        
        }
        { //::SireMM::CLJ14Group::usingGeometricCombiningRules
        
            typedef bool ( ::SireMM::CLJ14Group::*usingGeometricCombiningRules_function_type)(  ) const;
            usingGeometricCombiningRules_function_type usingGeometricCombiningRules_function_value( &::SireMM::CLJ14Group::usingGeometricCombiningRules );
            
            CLJ14Group_exposer.def( 
                "usingGeometricCombiningRules"
                , usingGeometricCombiningRules_function_value
                , bp::release_gil_policy()
                , "Return whether or not geometric combining rules are used" );
        
        }
        { //::SireMM::CLJ14Group::what
        
            typedef char const * ( ::SireMM::CLJ14Group::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::CLJ14Group::what );
            
            CLJ14Group_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJ14Group::wouldChangeProperties
        
            typedef bool ( ::SireMM::CLJ14Group::*wouldChangeProperties_function_type)( ::SireBase::PropertyMap const & ) const;
            wouldChangeProperties_function_type wouldChangeProperties_function_value( &::SireMM::CLJ14Group::wouldChangeProperties );
            
            CLJ14Group_exposer.def( 
                "wouldChangeProperties"
                , wouldChangeProperties_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "Return whether or not the passed property map would change properties that\nare used by this calculation" );
        
        }
        CLJ14Group_exposer.staticmethod( "typeName" );
        CLJ14Group_exposer.def( "__copy__", &__copy__);
        CLJ14Group_exposer.def( "__deepcopy__", &__copy__);
        CLJ14Group_exposer.def( "clone", &__copy__);
        CLJ14Group_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::CLJ14Group >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJ14Group_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::CLJ14Group >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJ14Group_exposer.def_pickle(sire_pickle_suite< ::SireMM::CLJ14Group >());
        CLJ14Group_exposer.def( "__str__", &__str__< ::SireMM::CLJ14Group > );
        CLJ14Group_exposer.def( "__repr__", &__str__< ::SireMM::CLJ14Group > );
    }

}
