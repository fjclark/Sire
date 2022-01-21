// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "InterCLJFFBase.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intercljff.h"

#include "intercljff.h"

SireFF::Inter2BFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> > __copy__(const SireFF::Inter2BFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> > &other){ return SireFF::Inter2BFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_InterCLJFFBase_class(){

    { //::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >
        typedef bp::class_< SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >, bp::bases< SireMM::CLJPotentialInterface<SireMM::InterCLJPotential>, SireFF::G1FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > InterCLJFFBase_exposer_t;
        InterCLJFFBase_exposer_t InterCLJFFBase_exposer = InterCLJFFBase_exposer_t( "InterCLJFFBase", "", bp::init< >("") );
        bp::scope InterCLJFFBase_scope( InterCLJFFBase_exposer );
        InterCLJFFBase_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        InterCLJFFBase_exposer.def( bp::init< SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > const & >(( bp::arg("other") ), "") );
        { //::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::components
        
            typedef SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef ::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::Components const & ( ::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*components_function_type)(  ) const;
            components_function_type components_function_value( &::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::components );
            
            InterCLJFFBase_exposer.def( 
                "components"
                , components_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::containsProperty
        
            typedef SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef bool ( ::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*containsProperty_function_type)( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::containsProperty );
            
            InterCLJFFBase_exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") )
                , "" );
        
        }
        { //::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::mustNowRecalculateFromScratch
        
            typedef SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*mustNowRecalculateFromScratch_function_type)(  ) ;
            mustNowRecalculateFromScratch_function_type mustNowRecalculateFromScratch_function_value( &::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::mustNowRecalculateFromScratch );
            
            InterCLJFFBase_exposer.def( 
                "mustNowRecalculateFromScratch"
                , mustNowRecalculateFromScratch_function_value
                , "" );
        
        }
        InterCLJFFBase_exposer.def( bp::self != bp::self );
        { //::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::operator=
        
            typedef SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef ::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > & ( ::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*assign_function_type)( ::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > const & ) ;
            assign_function_type assign_function_value( &::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::operator= );
            
            InterCLJFFBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        InterCLJFFBase_exposer.def( bp::self == bp::self );
        { //::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::properties
        
            typedef SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef ::SireBase::Properties const & ( ::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::properties );
            
            InterCLJFFBase_exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::property
        
            typedef SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef ::SireBase::Property const & ( ::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*property_function_type)( ::QString const & ) const;
            property_function_type property_function_value( &::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::property );
            
            InterCLJFFBase_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::setProperty
        
            typedef SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef bool ( ::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::setProperty );
            
            InterCLJFFBase_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("property") )
                , "" );
        
        }
        { //::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::typeName
        
            typedef SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::typeName );
            
            InterCLJFFBase_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::what
        
            typedef SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef char const * ( ::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::Inter2BFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::what );
            
            InterCLJFFBase_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        InterCLJFFBase_exposer.staticmethod( "typeName" );
        InterCLJFFBase_exposer.def( "__copy__", &__copy__);
        InterCLJFFBase_exposer.def( "__deepcopy__", &__copy__);
        InterCLJFFBase_exposer.def( "clone", &__copy__);
        InterCLJFFBase_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::Inter2BFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        InterCLJFFBase_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::Inter2BFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        InterCLJFFBase_exposer.def( "__getstate_manages_dict__", true);
        InterCLJFFBase_exposer.def( "__safe_for_unpickling__", true);
        InterCLJFFBase_exposer.def( "__setstate__", &__setstate__base64< ::SireFF::Inter2BFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> > > );
        InterCLJFFBase_exposer.def( "__getstate__", &__getstate__base64< ::SireFF::Inter2BFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> > > );
        InterCLJFFBase_exposer.def( "__str__", &__str__< ::SireFF::Inter2BFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> > > );
        InterCLJFFBase_exposer.def( "__repr__", &__str__< ::SireFF::Inter2BFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> > > );
        InterCLJFFBase_exposer.def( "__len__", &__len_count< ::SireFF::Inter2BFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> > > );
    }

}
