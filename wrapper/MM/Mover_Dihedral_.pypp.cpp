// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Mover_Dihedral_.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireCAS/expression.h"

#include "SireCAS/symbol.h"

#include "SireCAS/values.h"

#include "SireMol/errors.h"

#include "SireMol/molecule.h"

#include "SireMol/mover.hpp"

#include "SireMol/selector.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "dihedral.h"

#include "fouratomfunctions.h"

#include "selectordihedral.h"

#include <QDebug>

#include "dihedral.h"

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

SireMol::Mover<SireMM::Dihedral> __copy__(const SireMol::Mover<SireMM::Dihedral> &other){ return SireMol::Mover<SireMM::Dihedral>(other); }

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Mover_Dihedral__class(){

    { //::SireMol::Mover< SireMM::Dihedral >
        typedef bp::class_< SireMol::Mover< SireMM::Dihedral >, bp::bases< SireMol::MoverBase, SireMM::Dihedral, SireMol::MoleculeView, SireBase::Property > > Mover_Dihedral__exposer_t;
        Mover_Dihedral__exposer_t Mover_Dihedral__exposer = Mover_Dihedral__exposer_t( "Mover_Dihedral_", "", bp::init< >("") );
        bp::scope Mover_Dihedral__scope( Mover_Dihedral__exposer );
        Mover_Dihedral__exposer.def( bp::init< SireMM::Dihedral const & >(( bp::arg("view") ), "") );
        Mover_Dihedral__exposer.def( bp::init< SireMM::Dihedral const &, SireMol::AtomSelection const & >(( bp::arg("view"), bp::arg("movable_atoms") ), "") );
        Mover_Dihedral__exposer.def( bp::init< SireMol::Mover< SireMM::Dihedral > const & >(( bp::arg("other") ), "") );
        { //::SireMol::Mover< SireMM::Dihedral >::align
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*align_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) ;
            align_function_type align_function_value( &::SireMol::Mover< SireMM::Dihedral >::align );
            
            Mover_Dihedral__exposer.def( 
                "align"
                , align_function_value
                , ( bp::arg("other"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::align
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*align_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const &,::SireBase::PropertyMap const & ) ;
            align_function_type align_function_value( &::SireMol::Mover< SireMM::Dihedral >::align );
            
            Mover_Dihedral__exposer.def( 
                "align"
                , align_function_value
                , ( bp::arg("other"), bp::arg("map0"), bp::arg("map1") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::align
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*align_function_type)( ::SireMol::MoleculeView const &,::SireMol::AtomMatcher const &,::SireBase::PropertyMap const & ) ;
            align_function_type align_function_value( &::SireMol::Mover< SireMM::Dihedral >::align );
            
            Mover_Dihedral__exposer.def( 
                "align"
                , align_function_value
                , ( bp::arg("other"), bp::arg("matcher"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::align
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*align_function_type)( ::SireMol::MoleculeView const &,::SireMol::AtomMatcher const &,::SireBase::PropertyMap const &,::SireBase::PropertyMap const & ) ;
            align_function_type align_function_value( &::SireMol::Mover< SireMM::Dihedral >::align );
            
            Mover_Dihedral__exposer.def( 
                "align"
                , align_function_value
                , ( bp::arg("other"), bp::arg("matcher"), bp::arg("map0"), bp::arg("map1") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::change
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*change_function_type)( ::SireMol::BondID const &,::SireUnits::Dimension::Length,::SireBase::PropertyMap const & ) ;
            change_function_type change_function_value( &::SireMol::Mover< SireMM::Dihedral >::change );
            
            Mover_Dihedral__exposer.def( 
                "change"
                , change_function_value
                , ( bp::arg("bond"), bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::change
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*change_function_type)( ::SireMol::AngleID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            change_function_type change_function_value( &::SireMol::Mover< SireMM::Dihedral >::change );
            
            Mover_Dihedral__exposer.def( 
                "change"
                , change_function_value
                , ( bp::arg("angle"), bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::change
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*change_function_type)( ::SireMol::DihedralID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            change_function_type change_function_value( &::SireMol::Mover< SireMM::Dihedral >::change );
            
            Mover_Dihedral__exposer.def( 
                "change"
                , change_function_value
                , ( bp::arg("dihedral"), bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::change
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*change_function_type)( ::SireMol::BondID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            change_function_type change_function_value( &::SireMol::Mover< SireMM::Dihedral >::change );
            
            Mover_Dihedral__exposer.def( 
                "change"
                , change_function_value
                , ( bp::arg("bond"), bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::change
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*change_function_type)( ::SireMol::ImproperID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            change_function_type change_function_value( &::SireMol::Mover< SireMM::Dihedral >::change );
            
            Mover_Dihedral__exposer.def( 
                "change"
                , change_function_value
                , ( bp::arg("improper"), bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::changeFrame
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*changeFrame_function_type)( ::SireMaths::AxisSet const &,::SireMaths::AxisSet const &,::SireBase::PropertyMap const & ) ;
            changeFrame_function_type changeFrame_function_value( &::SireMol::Mover< SireMM::Dihedral >::changeFrame );
            
            Mover_Dihedral__exposer.def( 
                "changeFrame"
                , changeFrame_function_value
                , ( bp::arg("from_frame"), bp::arg("to_frame"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::commit
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMM::Dihedral ( ::SireMol::Mover< SireMM::Dihedral >::*commit_function_type)(  ) const;
            commit_function_type commit_function_value( &::SireMol::Mover< SireMM::Dihedral >::commit );
            
            Mover_Dihedral__exposer.def( 
                "commit"
                , commit_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::mapInto
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*mapInto_function_type)( ::SireMaths::AxisSet const &,::SireBase::PropertyMap const & ) ;
            mapInto_function_type mapInto_function_value( &::SireMol::Mover< SireMM::Dihedral >::mapInto );
            
            Mover_Dihedral__exposer.def( 
                "mapInto"
                , mapInto_function_value
                , ( bp::arg("axes"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::operator=
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*assign_function_type)( ::SireMol::Mover< SireMM::Dihedral > const & ) ;
            assign_function_type assign_function_value( &::SireMol::Mover< SireMM::Dihedral >::operator= );
            
            Mover_Dihedral__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::operator=
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*assign_function_type)( ::SireMM::Dihedral const & ) ;
            assign_function_type assign_function_value( &::SireMol::Mover< SireMM::Dihedral >::operator= );
            
            Mover_Dihedral__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::rotate
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*rotate_function_type)( ::SireMaths::Quaternion const &,::SireMaths::Vector const &,::SireBase::PropertyMap const & ) ;
            rotate_function_type rotate_function_value( &::SireMol::Mover< SireMM::Dihedral >::rotate );
            
            Mover_Dihedral__exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("quat"), bp::arg("point"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::rotate
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*rotate_function_type)( ::SireMaths::Matrix const &,::SireMaths::Vector const &,::SireBase::PropertyMap const & ) ;
            rotate_function_type rotate_function_value( &::SireMol::Mover< SireMM::Dihedral >::rotate );
            
            Mover_Dihedral__exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("rotmat"), bp::arg("point"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::set
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*set_function_type)( ::SireMol::BondID const &,::SireUnits::Dimension::Length,::SireBase::PropertyMap const & ) ;
            set_function_type set_function_value( &::SireMol::Mover< SireMM::Dihedral >::set );
            
            Mover_Dihedral__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("bond"), bp::arg("value"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::set
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*set_function_type)( ::SireMol::AngleID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            set_function_type set_function_value( &::SireMol::Mover< SireMM::Dihedral >::set );
            
            Mover_Dihedral__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("angle"), bp::arg("value"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::set
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*set_function_type)( ::SireMol::DihedralID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            set_function_type set_function_value( &::SireMol::Mover< SireMM::Dihedral >::set );
            
            Mover_Dihedral__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("dihedral"), bp::arg("value"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::set
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*set_function_type)( ::SireMol::ImproperID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            set_function_type set_function_value( &::SireMol::Mover< SireMM::Dihedral >::set );
            
            Mover_Dihedral__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("improper"), bp::arg("value"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::setAll
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*setAll_function_type)( ::SireMol::DihedralID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            setAll_function_type setAll_function_value( &::SireMol::Mover< SireMM::Dihedral >::setAll );
            
            Mover_Dihedral__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("dihedral"), bp::arg("value"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::toString
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::QString ( ::SireMol::Mover< SireMM::Dihedral >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::Mover< SireMM::Dihedral >::toString );
            
            Mover_Dihedral__exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::transform
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*transform_function_type)( ::SireMaths::Transform const &,::SireBase::PropertyMap const & ) ;
            transform_function_type transform_function_value( &::SireMol::Mover< SireMM::Dihedral >::transform );
            
            Mover_Dihedral__exposer.def( 
                "transform"
                , transform_function_value
                , ( bp::arg("transform"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::translate
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef ::SireMol::Mover< SireMM::Dihedral > & ( ::SireMol::Mover< SireMM::Dihedral >::*translate_function_type)( ::SireMaths::Vector const &,::SireBase::PropertyMap const & ) ;
            translate_function_type translate_function_value( &::SireMol::Mover< SireMM::Dihedral >::translate );
            
            Mover_Dihedral__exposer.def( 
                "translate"
                , translate_function_value
                , ( bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMM::Dihedral >::typeName
        
            typedef SireMol::Mover< SireMM::Dihedral > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Mover< SireMM::Dihedral >::typeName );
            
            Mover_Dihedral__exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Mover_Dihedral__exposer.staticmethod( "typeName" );
        Mover_Dihedral__exposer.def( "__copy__", &__copy__);
        Mover_Dihedral__exposer.def( "__deepcopy__", &__copy__);
        Mover_Dihedral__exposer.def( "clone", &__copy__);
        Mover_Dihedral__exposer.def( "__str__", &__str__< ::SireMol::Mover<SireMM::Dihedral> > );
        Mover_Dihedral__exposer.def( "__repr__", &__str__< ::SireMol::Mover<SireMM::Dihedral> > );
        Mover_Dihedral__exposer.def( "__len__", &__len_size< ::SireMol::Mover<SireMM::Dihedral> > );
    }

}
