// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Restraint3D.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/errors.h"

#include "SireCAS/expression.h"

#include "SireCAS/symbols.h"

#include "SireCAS/values.h"

#include "SireError/errors.h"

#include "SireFF/forcetable.h"

#include "SireMol/moleculedata.h"

#include "SireMol/molecules.h"

#include "SireMol/molid.h"

#include "SireMol/molnum.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "restraint.h"

#include "restraint.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Restraint3D_class(){

    { //::SireMM::Restraint3D
        typedef bp::class_< SireMM::Restraint3D, bp::bases< SireMM::Restraint, SireBase::Property >, boost::noncopyable > Restraint3D_exposer_t;
        Restraint3D_exposer_t Restraint3D_exposer = Restraint3D_exposer_t( "Restraint3D", "This is the base class of all restraints that operate in 3 dimensions,\nand so can thus return the force on the molecule caused by the restraint\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope Restraint3D_scope( Restraint3D_exposer );
        { //::SireMM::Restraint3D::force
        
            typedef void ( ::SireMM::Restraint3D::*force_function_type)( ::SireFF::MolForceTable &,double ) const;
            force_function_type force_function_value( &::SireMM::Restraint3D::force );
            
            Restraint3D_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("scale_force")=1 )
                , "" );
        
        }
        { //::SireMM::Restraint3D::force
        
            typedef void ( ::SireMM::Restraint3D::*force_function_type)( ::SireFF::ForceTable &,double ) const;
            force_function_type force_function_value( &::SireMM::Restraint3D::force );
            
            Restraint3D_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("scale_force")=1 )
                , "" );
        
        }
        { //::SireMM::Restraint3D::setSpace
        
            typedef void ( ::SireMM::Restraint3D::*setSpace_function_type)( ::SireVol::Space const & ) ;
            setSpace_function_type setSpace_function_value( &::SireMM::Restraint3D::setSpace );
            
            Restraint3D_exposer.def( 
                "setSpace"
                , setSpace_function_value
                , ( bp::arg("space") )
                , "Set the 3D space in which this restraint operates" );
        
        }
        { //::SireMM::Restraint3D::space
        
            typedef ::SireVol::Space const & ( ::SireMM::Restraint3D::*space_function_type)(  ) const;
            space_function_type space_function_value( &::SireMM::Restraint3D::space );
            
            Restraint3D_exposer.def( 
                "space"
                , space_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the 3D space in which this restraint operates" );
        
        }
        { //::SireMM::Restraint3D::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::Restraint3D::typeName );
            
            Restraint3D_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMM::Restraint3D::usesMoleculesIn
        
            typedef bool ( ::SireMM::Restraint3D::*usesMoleculesIn_function_type)( ::SireFF::ForceTable const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireMM::Restraint3D::usesMoleculesIn );
            
            Restraint3D_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("forcetable") )
                , "" );
        
        }
        { //::SireMM::Restraint3D::usesMoleculesIn
        
            typedef bool ( ::SireMM::Restraint3D::*usesMoleculesIn_function_type)( ::SireMol::Molecules const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireMM::Restraint3D::usesMoleculesIn );
            
            Restraint3D_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molecules") )
                , "" );
        
        }
        Restraint3D_exposer.staticmethod( "typeName" );
        Restraint3D_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::Restraint3D >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Restraint3D_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::Restraint3D >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Restraint3D_exposer.def_pickle(sire_pickle_suite< ::SireMM::Restraint3D >());
        Restraint3D_exposer.def( "__str__", &__str__< ::SireMM::Restraint3D > );
        Restraint3D_exposer.def( "__repr__", &__str__< ::SireMM::Restraint3D > );
    }

}
