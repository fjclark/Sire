// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Mover_SelectorBond_.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireBase/slice.h"

#include "SireCAS/symbol.h"

#include "SireCAS/values.h"

#include "SireID/index.h"

#include "SireMol/molecule.h"

#include "SireMol/mover.hpp"

#include "SireMol/selector.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "selectorbond.h"

#include "twoatomfunctions.h"

#include <QDebug>

#include "selectorbond.h"

#include "SireMaths/align.h"

#include "SireMaths/axisset.h"

#include "SireMaths/matrix.h"

#include "SireMaths/quaternion.h"

#include "SireMaths/rotate.h"

#include "SireMaths/vectorproperty.h"

#include "SireMol/errors.h"

#include "SireUnits/units.h"

#include "SireVol/coordgroup.h"

#include "SireVol/space.h"

#include "angleid.h"

#include "atomcoords.h"

#include "atommatcher.h"

#include "atommatchers.h"

#include "bondid.h"

#include "connectivity.h"

#include "dihedralid.h"

#include "improperid.h"

#include "mover.h"

#include "tostring.h"

#include "weightfunction.h"

#include "mover.h"

SireMol::Mover<SireMM::SelectorBond> __copy__(const SireMol::Mover<SireMM::SelectorBond> &other){ return SireMol::Mover<SireMM::SelectorBond>(other); }

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Mover_SelectorBond__class(){

    { //::SireMol::Mover< SireMM::SelectorBond >
        typedef bp::class_< SireMol::Mover< SireMM::SelectorBond >, bp::bases< SireMol::MoverBase, SireMM::SelectorBond, SireMol::MoleculeView, SireBase::Property > > Mover_SelectorBond__exposer_t;
        Mover_SelectorBond__exposer_t Mover_SelectorBond__exposer = Mover_SelectorBond__exposer_t( "Mover_SelectorBond_", "", bp::init< >("") );
        bp::scope Mover_SelectorBond__scope( Mover_SelectorBond__exposer );
        Mover_SelectorBond__exposer.def( bp::init< SireMM::SelectorBond const & >(( bp::arg("view") ), "") );
        Mover_SelectorBond__exposer.def( bp::init< SireMM::SelectorBond const &, SireMol::AtomSelection const & >(( bp::arg("view"), bp::arg("movable_atoms") ), "") );
        Mover_SelectorBond__exposer.def( bp::init< SireMol::Mover< SireMM::SelectorBond > const & >(( bp::arg("other") ), "") );
        { //::SireMol::Mover< SireMM::SelectorBond >::align
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*align_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) ;
            align_function_type align_function_value( &::SireMol::Mover< SireMM::SelectorBond >::align );
            
            Mover_SelectorBond__exposer.def( 
                "align"
                , align_function_value
                , ( bp::arg("other"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::align
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*align_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const &,::SireBase::PropertyMap const & ) ;
            align_function_type align_function_value( &::SireMol::Mover< SireMM::SelectorBond >::align );
            
            Mover_SelectorBond__exposer.def( 
                "align"
                , align_function_value
                , ( bp::arg("other"), bp::arg("map0"), bp::arg("map1") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::align
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*align_function_type)( ::SireMol::MoleculeView const &,::SireMol::AtomMatcher const &,::SireBase::PropertyMap const & ) ;
            align_function_type align_function_value( &::SireMol::Mover< SireMM::SelectorBond >::align );
            
            Mover_SelectorBond__exposer.def( 
                "align"
                , align_function_value
                , ( bp::arg("other"), bp::arg("matcher"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::align
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*align_function_type)( ::SireMol::MoleculeView const &,::SireMol::AtomMatcher const &,::SireBase::PropertyMap const &,::SireBase::PropertyMap const & ) ;
            align_function_type align_function_value( &::SireMol::Mover< SireMM::SelectorBond >::align );
            
            Mover_SelectorBond__exposer.def( 
                "align"
                , align_function_value
                , ( bp::arg("other"), bp::arg("matcher"), bp::arg("map0"), bp::arg("map1") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::change
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*change_function_type)( ::SireMol::BondID const &,::SireUnits::Dimension::Length,::SireBase::PropertyMap const & ) ;
            change_function_type change_function_value( &::SireMol::Mover< SireMM::SelectorBond >::change );
            
            Mover_SelectorBond__exposer.def( 
                "change"
                , change_function_value
                , ( bp::arg("bond"), bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::change
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*change_function_type)( ::SireMol::AngleID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            change_function_type change_function_value( &::SireMol::Mover< SireMM::SelectorBond >::change );
            
            Mover_SelectorBond__exposer.def( 
                "change"
                , change_function_value
                , ( bp::arg("angle"), bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::change
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*change_function_type)( ::SireMol::DihedralID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            change_function_type change_function_value( &::SireMol::Mover< SireMM::SelectorBond >::change );
            
            Mover_SelectorBond__exposer.def( 
                "change"
                , change_function_value
                , ( bp::arg("dihedral"), bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::change
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*change_function_type)( ::SireMol::BondID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            change_function_type change_function_value( &::SireMol::Mover< SireMM::SelectorBond >::change );
            
            Mover_SelectorBond__exposer.def( 
                "change"
                , change_function_value
                , ( bp::arg("bond"), bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::change
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*change_function_type)( ::SireMol::ImproperID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            change_function_type change_function_value( &::SireMol::Mover< SireMM::SelectorBond >::change );
            
            Mover_SelectorBond__exposer.def( 
                "change"
                , change_function_value
                , ( bp::arg("improper"), bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::changeFrame
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*changeFrame_function_type)( ::SireMaths::AxisSet const &,::SireMaths::AxisSet const &,::SireBase::PropertyMap const & ) ;
            changeFrame_function_type changeFrame_function_value( &::SireMol::Mover< SireMM::SelectorBond >::changeFrame );
            
            Mover_SelectorBond__exposer.def( 
                "changeFrame"
                , changeFrame_function_value
                , ( bp::arg("from_frame"), bp::arg("to_frame"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::commit
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMM::SelectorBond ( ::SireMol::Mover< SireMM::SelectorBond >::*commit_function_type)(  ) const;
            commit_function_type commit_function_value( &::SireMol::Mover< SireMM::SelectorBond >::commit );
            
            Mover_SelectorBond__exposer.def( 
                "commit"
                , commit_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::mapInto
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*mapInto_function_type)( ::SireMaths::AxisSet const &,::SireBase::PropertyMap const & ) ;
            mapInto_function_type mapInto_function_value( &::SireMol::Mover< SireMM::SelectorBond >::mapInto );
            
            Mover_SelectorBond__exposer.def( 
                "mapInto"
                , mapInto_function_value
                , ( bp::arg("axes"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::operator=
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*assign_function_type)( ::SireMol::Mover< SireMM::SelectorBond > const & ) ;
            assign_function_type assign_function_value( &::SireMol::Mover< SireMM::SelectorBond >::operator= );
            
            Mover_SelectorBond__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::operator=
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*assign_function_type)( ::SireMM::SelectorBond const & ) ;
            assign_function_type assign_function_value( &::SireMol::Mover< SireMM::SelectorBond >::operator= );
            
            Mover_SelectorBond__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::rotate
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*rotate_function_type)( ::SireMaths::Quaternion const &,::SireMaths::Vector const &,::SireBase::PropertyMap const & ) ;
            rotate_function_type rotate_function_value( &::SireMol::Mover< SireMM::SelectorBond >::rotate );
            
            Mover_SelectorBond__exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("quat"), bp::arg("point"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::rotate
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*rotate_function_type)( ::SireMaths::Matrix const &,::SireMaths::Vector const &,::SireBase::PropertyMap const & ) ;
            rotate_function_type rotate_function_value( &::SireMol::Mover< SireMM::SelectorBond >::rotate );
            
            Mover_SelectorBond__exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("rotmat"), bp::arg("point"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::set
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*set_function_type)( ::SireMol::BondID const &,::SireUnits::Dimension::Length,::SireBase::PropertyMap const & ) ;
            set_function_type set_function_value( &::SireMol::Mover< SireMM::SelectorBond >::set );
            
            Mover_SelectorBond__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("bond"), bp::arg("value"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::set
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*set_function_type)( ::SireMol::AngleID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            set_function_type set_function_value( &::SireMol::Mover< SireMM::SelectorBond >::set );
            
            Mover_SelectorBond__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("angle"), bp::arg("value"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::set
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*set_function_type)( ::SireMol::DihedralID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            set_function_type set_function_value( &::SireMol::Mover< SireMM::SelectorBond >::set );
            
            Mover_SelectorBond__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("dihedral"), bp::arg("value"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::set
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*set_function_type)( ::SireMol::ImproperID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            set_function_type set_function_value( &::SireMol::Mover< SireMM::SelectorBond >::set );
            
            Mover_SelectorBond__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("improper"), bp::arg("value"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::setAll
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*setAll_function_type)( ::SireMol::DihedralID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            setAll_function_type setAll_function_value( &::SireMol::Mover< SireMM::SelectorBond >::setAll );
            
            Mover_SelectorBond__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("dihedral"), bp::arg("value"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::toString
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::QString ( ::SireMol::Mover< SireMM::SelectorBond >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::Mover< SireMM::SelectorBond >::toString );
            
            Mover_SelectorBond__exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::transform
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*transform_function_type)( ::SireMaths::Transform const &,::SireBase::PropertyMap const & ) ;
            transform_function_type transform_function_value( &::SireMol::Mover< SireMM::SelectorBond >::transform );
            
            Mover_SelectorBond__exposer.def( 
                "transform"
                , transform_function_value
                , ( bp::arg("transform"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::translate
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef ::SireMol::Mover< SireMM::SelectorBond > & ( ::SireMol::Mover< SireMM::SelectorBond >::*translate_function_type)( ::SireMaths::Vector const &,::SireBase::PropertyMap const & ) ;
            translate_function_type translate_function_value( &::SireMol::Mover< SireMM::SelectorBond >::translate );
            
            Mover_SelectorBond__exposer.def( 
                "translate"
                , translate_function_value
                , ( bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::SelectorBond >::typeName
        
            typedef SireMol::Mover< SireMM::SelectorBond > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Mover< SireMM::SelectorBond >::typeName );
            
            Mover_SelectorBond__exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Mover_SelectorBond__exposer.staticmethod( "typeName" );
        Mover_SelectorBond__exposer.def( "__copy__", &__copy__);
        Mover_SelectorBond__exposer.def( "__deepcopy__", &__copy__);
        Mover_SelectorBond__exposer.def( "clone", &__copy__);
        Mover_SelectorBond__exposer.def( "__str__", &__str__< ::SireMol::Mover<SireMM::SelectorBond> > );
        Mover_SelectorBond__exposer.def( "__repr__", &__str__< ::SireMol::Mover<SireMM::SelectorBond> > );
        Mover_SelectorBond__exposer.def( "__len__", &__len_size< ::SireMol::Mover<SireMM::SelectorBond> > );
    }

}