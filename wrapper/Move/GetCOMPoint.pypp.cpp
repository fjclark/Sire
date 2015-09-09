// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "GetCOMPoint.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/vector.h"

#include "SireMol/atom.h"

#include "SireMol/atomcoords.h"

#include "SireMol/evaluator.h"

#include "SireMol/moleculeview.h"

#include "SireMol/mover.hpp"

#include "SireMol/selector.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "getpoint.h"

#include <QDebug>

#include "getpoint.h"

SireMove::GetCOMPoint __copy__(const SireMove::GetCOMPoint &other){ return SireMove::GetCOMPoint(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_GetCOMPoint_class(){

    { //::SireMove::GetCOMPoint
        typedef bp::class_< SireMove::GetCOMPoint, bp::bases< SireMove::GetPoint, SireBase::Property > > GetCOMPoint_exposer_t;
        GetCOMPoint_exposer_t GetCOMPoint_exposer = GetCOMPoint_exposer_t( "GetCOMPoint", bp::init< >() );
        bp::scope GetCOMPoint_scope( GetCOMPoint_exposer );
        GetCOMPoint_exposer.def( bp::init< SireMol::AtomID const & >(( bp::arg("atomid") )) );
        GetCOMPoint_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const & >(( bp::arg("atomid0"), bp::arg("atomid1") )) );
        GetCOMPoint_exposer.def( bp::init< QList< SireMol::AtomIdentifier > const & >(( bp::arg("atomids") )) );
        GetCOMPoint_exposer.def( bp::init< SireMove::GetCOMPoint const & >(( bp::arg("other") )) );
        { //::SireMove::GetCOMPoint::atomID
        
            typedef ::SireMol::AtomID const & ( ::SireMove::GetCOMPoint::*atomID_function_type)(  ) const;
            atomID_function_type atomID_function_value( &::SireMove::GetCOMPoint::atomID );
            
            GetCOMPoint_exposer.def( 
                "atomID"
                , atomID_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMove::GetCOMPoint::getPoint
        
            typedef ::SireMaths::Vector ( ::SireMove::GetCOMPoint::*getPoint_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            getPoint_function_type getPoint_function_value( &::SireMove::GetCOMPoint::getPoint );
            
            GetCOMPoint_exposer.def( 
                "getPoint"
                , getPoint_function_value
                , ( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() ) );
        
        }
        GetCOMPoint_exposer.def( bp::self != bp::self );
        { //::SireMove::GetCOMPoint::operator=
        
            typedef ::SireMove::GetCOMPoint & ( ::SireMove::GetCOMPoint::*assign_function_type)( ::SireMove::GetCOMPoint const & ) ;
            assign_function_type assign_function_value( &::SireMove::GetCOMPoint::operator= );
            
            GetCOMPoint_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        GetCOMPoint_exposer.def( bp::self == bp::self );
        { //::SireMove::GetCOMPoint::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::GetCOMPoint::typeName );
            
            GetCOMPoint_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        GetCOMPoint_exposer.staticmethod( "typeName" );
        GetCOMPoint_exposer.def( "__copy__", &__copy__);
        GetCOMPoint_exposer.def( "__deepcopy__", &__copy__);
        GetCOMPoint_exposer.def( "clone", &__copy__);
        GetCOMPoint_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::GetCOMPoint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GetCOMPoint_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::GetCOMPoint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GetCOMPoint_exposer.def( "__str__", &__str__< ::SireMove::GetCOMPoint > );
        GetCOMPoint_exposer.def( "__repr__", &__str__< ::SireMove::GetCOMPoint > );
    }

}
