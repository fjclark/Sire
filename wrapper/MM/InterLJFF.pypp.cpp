// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "InterLJFF.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "interljff.h"

#include "interljff.h"

SireFF::Inter2B3DFF<SireMM::LJPotentialInterface<SireMM::InterLJPotential> > __copy__(const SireFF::Inter2B3DFF<SireMM::LJPotentialInterface<SireMM::InterLJPotential> > &other){ return SireFF::Inter2B3DFF<SireMM::LJPotentialInterface<SireMM::InterLJPotential> >(other); }

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_InterLJFF_class(){

    { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >
        typedef bp::class_< SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >, bp::bases< SireFF::FF3D, SireFF::Inter2BFF<SireMM::LJPotentialInterface<SireMM::InterLJPotential> >, SireMM::LJPotentialInterface<SireMM::InterLJPotential>, SireFF::G1FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > InterLJFF_exposer_t;
        InterLJFF_exposer_t InterLJFF_exposer = InterLJFF_exposer_t( "InterLJFF", "", bp::init< >("") );
        bp::scope InterLJFF_scope( InterLJFF_exposer );
        InterLJFF_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        InterLJFF_exposer.def( bp::init< SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > const & >(( bp::arg("other") ), "") );
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::energy
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*energy_function_type)(  ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::energy );
            
            InterLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::energy
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*energy_function_type)( ::SireCAS::Symbol const & ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::energy );
            
            InterLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("component") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::energy
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*energy_function_type)( ::SireFF::EnergyTable &,double ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::energy );
            
            InterLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("scale_energy")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::energy
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*energy_function_type)( ::SireFF::EnergyTable &,::SireCAS::Symbol const &,double ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::energy );
            
            InterLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("symbol"), bp::arg("scale_energy")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::field );
            
            InterLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,::SireCAS::Symbol const &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::field );
            
            InterLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::field );
            
            InterLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("probe"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::field );
            
            InterLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::force
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*force_function_type)( ::SireFF::ForceTable &,double ) ;
            force_function_type force_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::force );
            
            InterLJFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("scale_force")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::force
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*force_function_type)( ::SireFF::ForceTable &,::SireCAS::Symbol const &,double ) ;
            force_function_type force_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::force );
            
            InterLJFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("symbol"), bp::arg("scale_force")=1 )
                , "" );
        
        }
        InterLJFF_exposer.def( bp::self != bp::self );
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::operator=
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > & ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*assign_function_type)( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > const & ) ;
            assign_function_type assign_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::operator= );
            
            InterLJFF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        InterLJFF_exposer.def( bp::self == bp::self );
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::packCoordinates
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*packCoordinates_function_type)(  ) ;
            packCoordinates_function_type packCoordinates_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::packCoordinates );
            
            InterLJFF_exposer.def( 
                "packCoordinates"
                , packCoordinates_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::potential );
            
            InterLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::potential );
            
            InterLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::potential );
            
            InterLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("probe"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::potential );
            
            InterLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::typeName
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::typeName );
            
            InterLJFF_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::what
        
            typedef SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > > exported_class_t;
            typedef char const * ( ::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::Inter2B3DFF< SireMM::LJPotentialInterface< SireMM::InterLJPotential > >::what );
            
            InterLJFF_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        InterLJFF_exposer.staticmethod( "typeName" );
        InterLJFF_exposer.def( "__copy__", &__copy__);
        InterLJFF_exposer.def( "__deepcopy__", &__copy__);
        InterLJFF_exposer.def( "clone", &__copy__);
        InterLJFF_exposer.def( "__str__", &__str__< ::SireFF::Inter2B3DFF<SireMM::LJPotentialInterface<SireMM::InterLJPotential> > > );
        InterLJFF_exposer.def( "__repr__", &__str__< ::SireFF::Inter2B3DFF<SireMM::LJPotentialInterface<SireMM::InterLJPotential> > > );
        InterLJFF_exposer.def( "__len__", &__len_count< ::SireFF::Inter2B3DFF<SireMM::LJPotentialInterface<SireMM::InterLJPotential> > > );
    }

}
