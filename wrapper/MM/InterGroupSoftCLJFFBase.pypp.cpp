// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "InterGroupSoftCLJFFBase.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intersoftcljff.h"

#include "intersoftcljff.h"

SireFF::Inter2B2GFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > __copy__(const SireFF::Inter2B2GFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > &other){ return SireFF::Inter2B2GFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> >(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_InterGroupSoftCLJFFBase_class(){

    { //::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >
        typedef bp::class_< SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >, bp::bases< SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential>, SireMM::CLJPotentialInterface<SireMM::InterSoftCLJPotential>, SireFF::G2FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > InterGroupSoftCLJFFBase_exposer_t;
        InterGroupSoftCLJFFBase_exposer_t InterGroupSoftCLJFFBase_exposer = InterGroupSoftCLJFFBase_exposer_t( "InterGroupSoftCLJFFBase", "", bp::init< >("") );
        bp::scope InterGroupSoftCLJFFBase_scope( InterGroupSoftCLJFFBase_exposer );
        InterGroupSoftCLJFFBase_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        InterGroupSoftCLJFFBase_exposer.def( bp::init< SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > const & >(( bp::arg("other") ), "") );
        { //::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::components
        
            typedef SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef ::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::Components const & ( ::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*components_function_type)(  ) const;
            components_function_type components_function_value( &::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::components );
            
            InterGroupSoftCLJFFBase_exposer.def( 
                "components"
                , components_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::containsProperty
        
            typedef SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef bool ( ::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*containsProperty_function_type)( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::containsProperty );
            
            InterGroupSoftCLJFFBase_exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::mustNowRecalculateFromScratch
        
            typedef SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*mustNowRecalculateFromScratch_function_type)(  ) ;
            mustNowRecalculateFromScratch_function_type mustNowRecalculateFromScratch_function_value( &::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::mustNowRecalculateFromScratch );
            
            InterGroupSoftCLJFFBase_exposer.def( 
                "mustNowRecalculateFromScratch"
                , mustNowRecalculateFromScratch_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        InterGroupSoftCLJFFBase_exposer.def( bp::self != bp::self );
        { //::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::operator=
        
            typedef SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef ::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > & ( ::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*assign_function_type)( ::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > const & ) ;
            assign_function_type assign_function_value( &::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::operator= );
            
            InterGroupSoftCLJFFBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        InterGroupSoftCLJFFBase_exposer.def( bp::self == bp::self );
        { //::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::properties
        
            typedef SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef ::SireBase::Properties const & ( ::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::properties );
            
            InterGroupSoftCLJFFBase_exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::property
        
            typedef SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef ::SireBase::Property const & ( ::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*property_function_type)( ::QString const & ) const;
            property_function_type property_function_value( &::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::property );
            
            InterGroupSoftCLJFFBase_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::setProperty
        
            typedef SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef bool ( ::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::setProperty );
            
            InterGroupSoftCLJFFBase_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("property") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::typeName
        
            typedef SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::typeName );
            
            InterGroupSoftCLJFFBase_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::what
        
            typedef SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef char const * ( ::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::Inter2B2GFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::what );
            
            InterGroupSoftCLJFFBase_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        InterGroupSoftCLJFFBase_exposer.staticmethod( "typeName" );
        InterGroupSoftCLJFFBase_exposer.def( "__copy__", &__copy__);
        InterGroupSoftCLJFFBase_exposer.def( "__deepcopy__", &__copy__);
        InterGroupSoftCLJFFBase_exposer.def( "clone", &__copy__);
        InterGroupSoftCLJFFBase_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::Inter2B2GFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        InterGroupSoftCLJFFBase_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::Inter2B2GFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        InterGroupSoftCLJFFBase_exposer.def_pickle(sire_pickle_suite< ::SireFF::Inter2B2GFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > >());
        InterGroupSoftCLJFFBase_exposer.def( "__str__", &__str__< ::SireFF::Inter2B2GFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > > );
        InterGroupSoftCLJFFBase_exposer.def( "__repr__", &__str__< ::SireFF::Inter2B2GFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > > );
        InterGroupSoftCLJFFBase_exposer.def( "__len__", &__len_count< ::SireFF::Inter2B2GFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > > );
    }

}
