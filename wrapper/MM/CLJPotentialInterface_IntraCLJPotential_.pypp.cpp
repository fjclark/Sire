// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "CLJPotentialInterface_IntraCLJPotential_.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intracljff.h"

#include "intracljff.h"

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential>&){ return "SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential>";}

void register_CLJPotentialInterface_IntraCLJPotential__class(){

    { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >
        typedef bp::class_< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >, boost::noncopyable > CLJPotentialInterface_IntraCLJPotential__exposer_t;
        CLJPotentialInterface_IntraCLJPotential__exposer_t CLJPotentialInterface_IntraCLJPotential__exposer = CLJPotentialInterface_IntraCLJPotential__exposer_t( "CLJPotentialInterface_IntraCLJPotential_", "", bp::no_init );
        bp::scope CLJPotentialInterface_IntraCLJPotential__scope( CLJPotentialInterface_IntraCLJPotential__exposer );
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::combiningRules
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef ::QString const & ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*combiningRules_function_type)(  ) const;
            combiningRules_function_type combiningRules_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::combiningRules );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "combiningRules"
                , combiningRules_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::containsProperty
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*containsProperty_function_type)( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::containsProperty );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") )
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::parameters
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef ::SireMM::IntraCLJPotential::ParameterNames ( *parameters_function_type )(  );
            parameters_function_type parameters_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::parameters );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "parameters"
                , parameters_function_value
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::properties
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef ::SireBase::Properties const & ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::properties );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::property
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef ::SireBase::Property const & ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*property_function_type)( ::QString const & ) const;
            property_function_type property_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::property );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::reactionFieldDielectric
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef double ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*reactionFieldDielectric_function_type)(  ) const;
            reactionFieldDielectric_function_type reactionFieldDielectric_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::reactionFieldDielectric );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "reactionFieldDielectric"
                , reactionFieldDielectric_function_value
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setCombiningRules
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*setCombiningRules_function_type)( ::QString const & ) ;
            setCombiningRules_function_type setCombiningRules_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setCombiningRules );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "setCombiningRules"
                , setCombiningRules_function_value
                , ( bp::arg("combiningrules") )
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setProperty
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setProperty );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setReactionFieldDielectric
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*setReactionFieldDielectric_function_type)( double ) ;
            setReactionFieldDielectric_function_type setReactionFieldDielectric_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setReactionFieldDielectric );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "setReactionFieldDielectric"
                , setReactionFieldDielectric_function_value
                , ( bp::arg("dielectric") )
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setShiftElectrostatics
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*setShiftElectrostatics_function_type)( bool ) ;
            setShiftElectrostatics_function_type setShiftElectrostatics_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setShiftElectrostatics );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "setShiftElectrostatics"
                , setShiftElectrostatics_function_value
                , ( bp::arg("switchelectro") )
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setSpace
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*setSpace_function_type)( ::SireVol::Space const & ) ;
            setSpace_function_type setSpace_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setSpace );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "setSpace"
                , setSpace_function_value
                , ( bp::arg("new_space") )
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setSwitchingFunction
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*setSwitchingFunction_function_type)( ::SireMM::SwitchingFunction const & ) ;
            setSwitchingFunction_function_type setSwitchingFunction_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setSwitchingFunction );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "setSwitchingFunction"
                , setSwitchingFunction_function_value
                , ( bp::arg("new_switchfunc") )
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setUseAtomisticCutoff
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*setUseAtomisticCutoff_function_type)( bool ) ;
            setUseAtomisticCutoff_function_type setUseAtomisticCutoff_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setUseAtomisticCutoff );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "setUseAtomisticCutoff"
                , setUseAtomisticCutoff_function_value
                , ( bp::arg("switchatomistic") )
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setUseGroupCutoff
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*setUseGroupCutoff_function_type)( bool ) ;
            setUseGroupCutoff_function_type setUseGroupCutoff_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setUseGroupCutoff );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "setUseGroupCutoff"
                , setUseGroupCutoff_function_value
                , ( bp::arg("switchgroup") )
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setUseReactionField
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*setUseReactionField_function_type)( bool ) ;
            setUseReactionField_function_type setUseReactionField_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::setUseReactionField );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "setUseReactionField"
                , setUseReactionField_function_value
                , ( bp::arg("switchrf") )
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::shiftElectrostatics
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*shiftElectrostatics_function_type)(  ) const;
            shiftElectrostatics_function_type shiftElectrostatics_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::shiftElectrostatics );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "shiftElectrostatics"
                , shiftElectrostatics_function_value
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::space
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef ::SireVol::Space const & ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*space_function_type)(  ) const;
            space_function_type space_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::space );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "space"
                , space_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::switchingFunction
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef ::SireMM::SwitchingFunction const & ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*switchingFunction_function_type)(  ) const;
            switchingFunction_function_type switchingFunction_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::switchingFunction );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "switchingFunction"
                , switchingFunction_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::useAtomisticCutoff
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*useAtomisticCutoff_function_type)(  ) const;
            useAtomisticCutoff_function_type useAtomisticCutoff_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::useAtomisticCutoff );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "useAtomisticCutoff"
                , useAtomisticCutoff_function_value
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::useGroupCutoff
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*useGroupCutoff_function_type)(  ) const;
            useGroupCutoff_function_type useGroupCutoff_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::useGroupCutoff );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "useGroupCutoff"
                , useGroupCutoff_function_value
                , "" );
        
        }
        { //::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::useReactionField
        
            typedef SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::*useReactionField_function_type)(  ) const;
            useReactionField_function_type useReactionField_function_value( &::SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential >::useReactionField );
            
            CLJPotentialInterface_IntraCLJPotential__exposer.def( 
                "useReactionField"
                , useReactionField_function_value
                , "" );
        
        }
        CLJPotentialInterface_IntraCLJPotential__exposer.staticmethod( "parameters" );
        CLJPotentialInterface_IntraCLJPotential__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJPotentialInterface_IntraCLJPotential__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJPotentialInterface_IntraCLJPotential__exposer.def_pickle(sire_pickle_suite< ::SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> >());
        CLJPotentialInterface_IntraCLJPotential__exposer.def( "__str__", &pvt_get_name);
        CLJPotentialInterface_IntraCLJPotential__exposer.def( "__repr__", &pvt_get_name);
    }

}
