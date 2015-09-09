// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "IntraGroupCoulombFFBase.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intracoulombff.h"

#include "intracoulombff.h"

SireFF::Intra2B2GFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> > __copy__(const SireFF::Intra2B2GFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> > &other){ return SireFF::Intra2B2GFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> >(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_IntraGroupCoulombFFBase_class(){

    { //::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >
        typedef bp::class_< SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >, bp::bases< SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential>, SireFF::G2FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > IntraGroupCoulombFFBase_exposer_t;
        IntraGroupCoulombFFBase_exposer_t IntraGroupCoulombFFBase_exposer = IntraGroupCoulombFFBase_exposer_t( "IntraGroupCoulombFFBase", bp::init< >() );
        bp::scope IntraGroupCoulombFFBase_scope( IntraGroupCoulombFFBase_exposer );
        IntraGroupCoulombFFBase_exposer.def( bp::init< QString const & >(( bp::arg("name") )) );
        IntraGroupCoulombFFBase_exposer.def( bp::init< SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > > const & >(( bp::arg("other") )) );
        { //::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::components
        
            typedef SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > > exported_class_t;
            typedef ::SireMM::CoulombComponent const & ( ::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::*components_function_type)(  ) const;
            components_function_type components_function_value( &::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::components );
            
            IntraGroupCoulombFFBase_exposer.def( 
                "components"
                , components_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::containsProperty
        
            typedef SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > > exported_class_t;
            typedef bool ( ::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::*containsProperty_function_type)( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::containsProperty );
            
            IntraGroupCoulombFFBase_exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") ) );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::mustNowRecalculateFromScratch
        
            typedef SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::*mustNowRecalculateFromScratch_function_type)(  ) ;
            mustNowRecalculateFromScratch_function_type mustNowRecalculateFromScratch_function_value( &::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::mustNowRecalculateFromScratch );
            
            IntraGroupCoulombFFBase_exposer.def( 
                "mustNowRecalculateFromScratch"
                , mustNowRecalculateFromScratch_function_value );
        
        }
        IntraGroupCoulombFFBase_exposer.def( bp::self != bp::self );
        { //::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::operator=
        
            typedef SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > > exported_class_t;
            typedef ::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > > & ( ::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::*assign_function_type)( ::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > > const & ) ;
            assign_function_type assign_function_value( &::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::operator= );
            
            IntraGroupCoulombFFBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        IntraGroupCoulombFFBase_exposer.def( bp::self == bp::self );
        { //::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::properties
        
            typedef SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > > exported_class_t;
            typedef ::SireBase::Properties const & ( ::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::properties );
            
            IntraGroupCoulombFFBase_exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::property
        
            typedef SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > > exported_class_t;
            typedef ::SireBase::Property const & ( ::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::*property_function_type)( ::QString const & ) const;
            property_function_type property_function_value( &::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::property );
            
            IntraGroupCoulombFFBase_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::setProperty
        
            typedef SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > > exported_class_t;
            typedef bool ( ::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::setProperty );
            
            IntraGroupCoulombFFBase_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("property") ) );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::typeName
        
            typedef SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::typeName );
            
            IntraGroupCoulombFFBase_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::what
        
            typedef SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > > exported_class_t;
            typedef char const * ( ::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::Intra2B2GFF< SireMM::CoulombPotentialInterface< SireMM::IntraCoulombPotential > >::what );
            
            IntraGroupCoulombFFBase_exposer.def( 
                "what"
                , what_function_value );
        
        }
        IntraGroupCoulombFFBase_exposer.staticmethod( "typeName" );
        IntraGroupCoulombFFBase_exposer.def( "__copy__", &__copy__);
        IntraGroupCoulombFFBase_exposer.def( "__deepcopy__", &__copy__);
        IntraGroupCoulombFFBase_exposer.def( "clone", &__copy__);
        IntraGroupCoulombFFBase_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::Intra2B2GFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IntraGroupCoulombFFBase_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::Intra2B2GFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IntraGroupCoulombFFBase_exposer.def( "__str__", &__str__< ::SireFF::Intra2B2GFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> > > );
        IntraGroupCoulombFFBase_exposer.def( "__repr__", &__str__< ::SireFF::Intra2B2GFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> > > );
        IntraGroupCoulombFFBase_exposer.def( "__len__", &__len_count< ::SireFF::Intra2B2GFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> > > );
    }

}
