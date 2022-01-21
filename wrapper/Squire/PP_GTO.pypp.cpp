// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "PP_GTO.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/array2d.hpp"

#include "SireError/errors.h"

#include "SireMaths/boys.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "pgto.h"

#include "pointcharge.h"

#include "pointdipole.h"

#include "pgto.h"

Squire::PP_GTO __copy__(const Squire::PP_GTO &other){ return Squire::PP_GTO(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_PP_GTO_class(){

    { //::Squire::PP_GTO
        typedef bp::class_< Squire::PP_GTO, bp::bases< Squire::GTOPair, Squire::ShellPair, SireBase::Property > > PP_GTO_exposer_t;
        PP_GTO_exposer_t PP_GTO_exposer = PP_GTO_exposer_t( "PP_GTO", "This is a combined P-P GTO shell pair", bp::init< >("Constructor") );
        bp::scope PP_GTO_scope( PP_GTO_exposer );
        PP_GTO_exposer.def( bp::init< SireMaths::Vector const &, Squire::P_GTO const &, SireMaths::Vector const &, Squire::P_GTO const & >(( bp::arg("A"), bp::arg("a"), bp::arg("B"), bp::arg("b") ), "Construct combining orbital a at position A with orbital b at\nposition B") );
        PP_GTO_exposer.def( bp::init< Squire::PP_GTO const & >(( bp::arg("other") ), "Copy constructor") );
        { //::Squire::PP_GTO::P_minus_A
        
            typedef ::SireMaths::Vector const & ( ::Squire::PP_GTO::*P_minus_A_function_type)(  ) const;
            P_minus_A_function_type P_minus_A_function_value( &::Squire::PP_GTO::P_minus_A );
            
            PP_GTO_exposer.def( 
                "P_minus_A"
                , P_minus_A_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the vector from the center of the first P orbital shell\nto the center of mass of the combined gaussian" );
        
        }
        { //::Squire::PP_GTO::P_minus_B
        
            typedef ::SireMaths::Vector const & ( ::Squire::PP_GTO::*P_minus_B_function_type)(  ) const;
            P_minus_B_function_type P_minus_B_function_value( &::Squire::PP_GTO::P_minus_B );
            
            PP_GTO_exposer.def( 
                "P_minus_B"
                , P_minus_B_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the vector from the center of the second P orbital shell\nto the center of mass of the combined gaussian" );
        
        }
        { //::Squire::PP_GTO::Q_minus_C
        
            typedef ::SireMaths::Vector const & ( ::Squire::PP_GTO::*Q_minus_C_function_type)(  ) const;
            Q_minus_C_function_type Q_minus_C_function_value( &::Squire::PP_GTO::Q_minus_C );
            
            PP_GTO_exposer.def( 
                "Q_minus_C"
                , Q_minus_C_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Synonym for P_minus_A" );
        
        }
        { //::Squire::PP_GTO::Q_minus_D
        
            typedef ::SireMaths::Vector const & ( ::Squire::PP_GTO::*Q_minus_D_function_type)(  ) const;
            Q_minus_D_function_type Q_minus_D_function_value( &::Squire::PP_GTO::Q_minus_D );
            
            PP_GTO_exposer.def( 
                "Q_minus_D"
                , Q_minus_D_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Synonym for Q_minus_D" );
        
        }
        { //::Squire::PP_GTO::angularMomentum0
        
            typedef int ( ::Squire::PP_GTO::*angularMomentum0_function_type)(  ) const;
            angularMomentum0_function_type angularMomentum0_function_value( &::Squire::PP_GTO::angularMomentum0 );
            
            PP_GTO_exposer.def( 
                "angularMomentum0"
                , angularMomentum0_function_value
                , "Return the angular momentum of the first GTO shell in this pair" );
        
        }
        { //::Squire::PP_GTO::angularMomentum1
        
            typedef int ( ::Squire::PP_GTO::*angularMomentum1_function_type)(  ) const;
            angularMomentum1_function_type angularMomentum1_function_value( &::Squire::PP_GTO::angularMomentum1 );
            
            PP_GTO_exposer.def( 
                "angularMomentum1"
                , angularMomentum1_function_value
                , "Return the angular momentum of the second GTO shell in this pair" );
        
        }
        { //::Squire::PP_GTO::nOrbitals0
        
            typedef int ( ::Squire::PP_GTO::*nOrbitals0_function_type)(  ) const;
            nOrbitals0_function_type nOrbitals0_function_value( &::Squire::PP_GTO::nOrbitals0 );
            
            PP_GTO_exposer.def( 
                "nOrbitals0"
                , nOrbitals0_function_value
                , "Return the number of orbitals in the first GTO shell in this pair" );
        
        }
        { //::Squire::PP_GTO::nOrbitals1
        
            typedef int ( ::Squire::PP_GTO::*nOrbitals1_function_type)(  ) const;
            nOrbitals1_function_type nOrbitals1_function_value( &::Squire::PP_GTO::nOrbitals1 );
            
            PP_GTO_exposer.def( 
                "nOrbitals1"
                , nOrbitals1_function_value
                , "Return the number of orbitals in the second GTO shell in this pair" );
        
        }
        PP_GTO_exposer.def( bp::self != bp::self );
        { //::Squire::PP_GTO::operator=
        
            typedef ::Squire::PP_GTO & ( ::Squire::PP_GTO::*assign_function_type)( ::Squire::PP_GTO const & ) ;
            assign_function_type assign_function_value( &::Squire::PP_GTO::operator= );
            
            PP_GTO_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        PP_GTO_exposer.def( bp::self == bp::self );
        { //::Squire::PP_GTO::scale
        
            typedef double ( ::Squire::PP_GTO::*scale_function_type)(  ) const;
            scale_function_type scale_function_value( &::Squire::PP_GTO::scale );
            
            PP_GTO_exposer.def( 
                "scale"
                , scale_function_value
                , "Return the additional scaling constant needed to normalise the\nintegrals involving this shell-pair" );
        
        }
        { //::Squire::PP_GTO::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Squire::PP_GTO::typeName );
            
            PP_GTO_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        PP_GTO_exposer.staticmethod( "typeName" );
        PP_GTO_exposer.def( "__copy__", &__copy__);
        PP_GTO_exposer.def( "__deepcopy__", &__copy__);
        PP_GTO_exposer.def( "clone", &__copy__);
        PP_GTO_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Squire::PP_GTO >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PP_GTO_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Squire::PP_GTO >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PP_GTO_exposer.def( "__getstate_manages_dict__", true);
        PP_GTO_exposer.def( "__safe_for_unpickling__", true);
        PP_GTO_exposer.def( "__setstate__", &__setstate__base64< ::Squire::PP_GTO > );
        PP_GTO_exposer.def( "__getstate__", &__getstate__base64< ::Squire::PP_GTO > );
        PP_GTO_exposer.def( "__str__", &__str__< ::Squire::PP_GTO > );
        PP_GTO_exposer.def( "__repr__", &__str__< ::Squire::PP_GTO > );
    }

}
