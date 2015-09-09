// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "QMChargeConstraint.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/atomcharges.h"

#include "SireMol/mgname.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/molecules.h"

#include "SireMol/moleditor.h"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/delta.h"

#include "SireSystem/system.h"

#include "qmchargeconstraint.h"

#include "qmchargeconstraint.h"

Squire::QMChargeConstraint __copy__(const Squire::QMChargeConstraint &other){ return Squire::QMChargeConstraint(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_QMChargeConstraint_class(){

    { //::Squire::QMChargeConstraint
        typedef bp::class_< Squire::QMChargeConstraint, bp::bases< SireSystem::ChargeConstraint, SireSystem::MoleculeConstraint, SireSystem::Constraint, SireBase::Property > > QMChargeConstraint_exposer_t;
        QMChargeConstraint_exposer_t QMChargeConstraint_exposer = QMChargeConstraint_exposer_t( "QMChargeConstraint", bp::init< >() );
        bp::scope QMChargeConstraint_scope( QMChargeConstraint_exposer );
        QMChargeConstraint_exposer.def( bp::init< SireMol::MoleculeGroup const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() )) );
        QMChargeConstraint_exposer.def( bp::init< SireMol::MoleculeGroup const &, Squire::QMChargeCalculator const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("chargecalculator"), bp::arg("map")=SireBase::PropertyMap() )) );
        QMChargeConstraint_exposer.def( bp::init< Squire::QMChargeConstraint const & >(( bp::arg("other") )) );
        { //::Squire::QMChargeConstraint::chargeCalculator
        
            typedef ::Squire::QMChargeCalculator const & ( ::Squire::QMChargeConstraint::*chargeCalculator_function_type)(  ) const;
            chargeCalculator_function_type chargeCalculator_function_value( &::Squire::QMChargeConstraint::chargeCalculator );
            
            QMChargeConstraint_exposer.def( 
                "chargeCalculator"
                , chargeCalculator_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        QMChargeConstraint_exposer.def( bp::self != bp::self );
        { //::Squire::QMChargeConstraint::operator=
        
            typedef ::Squire::QMChargeConstraint & ( ::Squire::QMChargeConstraint::*assign_function_type)( ::Squire::QMChargeConstraint const & ) ;
            assign_function_type assign_function_value( &::Squire::QMChargeConstraint::operator= );
            
            QMChargeConstraint_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        QMChargeConstraint_exposer.def( bp::self == bp::self );
        { //::Squire::QMChargeConstraint::setChargeCalculator
        
            typedef void ( ::Squire::QMChargeConstraint::*setChargeCalculator_function_type)( ::Squire::QMChargeCalculator const & ) ;
            setChargeCalculator_function_type setChargeCalculator_function_value( &::Squire::QMChargeConstraint::setChargeCalculator );
            
            QMChargeConstraint_exposer.def( 
                "setChargeCalculator"
                , setChargeCalculator_function_value
                , ( bp::arg("chargecalculator") ) );
        
        }
        { //::Squire::QMChargeConstraint::toString
        
            typedef ::QString ( ::Squire::QMChargeConstraint::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::Squire::QMChargeConstraint::toString );
            
            QMChargeConstraint_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::Squire::QMChargeConstraint::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Squire::QMChargeConstraint::typeName );
            
            QMChargeConstraint_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        QMChargeConstraint_exposer.staticmethod( "typeName" );
        QMChargeConstraint_exposer.def( "__copy__", &__copy__);
        QMChargeConstraint_exposer.def( "__deepcopy__", &__copy__);
        QMChargeConstraint_exposer.def( "clone", &__copy__);
        QMChargeConstraint_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Squire::QMChargeConstraint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        QMChargeConstraint_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Squire::QMChargeConstraint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        QMChargeConstraint_exposer.def( "__str__", &__str__< ::Squire::QMChargeConstraint > );
        QMChargeConstraint_exposer.def( "__repr__", &__str__< ::Squire::QMChargeConstraint > );
    }

}
