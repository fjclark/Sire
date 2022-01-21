// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "IntraLJFFBase.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intraljff.h"

#include "intraljff.h"

SireFF::Intra2BFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > __copy__(const SireFF::Intra2BFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > &other){ return SireFF::Intra2BFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> >(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_IntraLJFFBase_class(){

    { //::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >
        typedef bp::class_< SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >, bp::bases< SireMM::LJPotentialInterface<SireMM::IntraLJPotential>, SireFF::G1FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > IntraLJFFBase_exposer_t;
        IntraLJFFBase_exposer_t IntraLJFFBase_exposer = IntraLJFFBase_exposer_t( "IntraLJFFBase", "", bp::init< >("") );
        bp::scope IntraLJFFBase_scope( IntraLJFFBase_exposer );
        IntraLJFFBase_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        IntraLJFFBase_exposer.def( bp::init< SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > const & >(( bp::arg("other") ), "") );
        { //::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::components
        
            typedef SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef ::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::Components const & ( ::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*components_function_type)(  ) const;
            components_function_type components_function_value( &::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::components );
            
            IntraLJFFBase_exposer.def( 
                "components"
                , components_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::containsProperty
        
            typedef SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef bool ( ::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*containsProperty_function_type)( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::containsProperty );
            
            IntraLJFFBase_exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") )
                , "" );
        
        }
        { //::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::mustNowRecalculateFromScratch
        
            typedef SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*mustNowRecalculateFromScratch_function_type)(  ) ;
            mustNowRecalculateFromScratch_function_type mustNowRecalculateFromScratch_function_value( &::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::mustNowRecalculateFromScratch );
            
            IntraLJFFBase_exposer.def( 
                "mustNowRecalculateFromScratch"
                , mustNowRecalculateFromScratch_function_value
                , "" );
        
        }
        IntraLJFFBase_exposer.def( bp::self != bp::self );
        { //::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::operator=
        
            typedef SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef ::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > & ( ::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*assign_function_type)( ::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > const & ) ;
            assign_function_type assign_function_value( &::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::operator= );
            
            IntraLJFFBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        IntraLJFFBase_exposer.def( bp::self == bp::self );
        { //::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::properties
        
            typedef SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef ::SireBase::Properties const & ( ::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::properties );
            
            IntraLJFFBase_exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::property
        
            typedef SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef ::SireBase::Property const & ( ::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*property_function_type)( ::QString const & ) const;
            property_function_type property_function_value( &::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::property );
            
            IntraLJFFBase_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::setProperty
        
            typedef SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef bool ( ::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::setProperty );
            
            IntraLJFFBase_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("property") )
                , "" );
        
        }
        { //::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::typeName
        
            typedef SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::typeName );
            
            IntraLJFFBase_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::what
        
            typedef SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > > exported_class_t;
            typedef char const * ( ::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::Intra2BFF< SireMM::LJPotentialInterface< SireMM::IntraLJPotential > >::what );
            
            IntraLJFFBase_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        IntraLJFFBase_exposer.staticmethod( "typeName" );
        IntraLJFFBase_exposer.def( "__copy__", &__copy__);
        IntraLJFFBase_exposer.def( "__deepcopy__", &__copy__);
        IntraLJFFBase_exposer.def( "clone", &__copy__);
        IntraLJFFBase_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::Intra2BFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IntraLJFFBase_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::Intra2BFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IntraLJFFBase_exposer.def( "__setstate__", &__setstate__base64< ::SireFF::Intra2BFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > > );
        IntraLJFFBase_exposer.def( "__getstate__", &__getstate__base64< ::SireFF::Intra2BFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > > );
        IntraLJFFBase_exposer.def( "__str__", &__str__< ::SireFF::Intra2BFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > > );
        IntraLJFFBase_exposer.def( "__repr__", &__str__< ::SireFF::Intra2BFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > > );
        IntraLJFFBase_exposer.def( "__len__", &__len_count< ::SireFF::Intra2BFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> > > );
    }

}
