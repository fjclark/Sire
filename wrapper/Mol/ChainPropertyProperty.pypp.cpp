// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ChainPropertyProperty.pypp.hpp"

namespace bp = boost::python;

#include "chainproperty.hpp"

#include "chainproperty.hpp"

#include "SireMaths/vector.h"

#include "SireMol/moleculeview.h"

SireMol::ChainProperty<SireBase::PropPtr<SireBase::Property> > __copy__(const SireMol::ChainProperty<SireBase::PropPtr<SireBase::Property> > &other){ return SireMol::ChainProperty<SireBase::PropPtr<SireBase::Property> >(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_ChainPropertyProperty_class(){

    { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >
        typedef bp::class_< SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >, bp::bases< SireMol::ChainProp, SireMol::MolViewProperty, SireBase::Property > > ChainPropertyProperty_exposer_t;
        ChainPropertyProperty_exposer_t ChainPropertyProperty_exposer = ChainPropertyProperty_exposer_t( "ChainPropertyProperty", "", bp::init< >("") );
        bp::scope ChainPropertyProperty_scope( ChainPropertyProperty_exposer );
        ChainPropertyProperty_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") ), "") );
        ChainPropertyProperty_exposer.def( bp::init< QVector< SireBase::PropPtr< SireBase::Property > > const & >(( bp::arg("values") ), "") );
        ChainPropertyProperty_exposer.def( bp::init< SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > const & >(( bp::arg("other") ), "") );
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::array
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::QVector< SireBase::PropPtr< SireBase::Property > > const & ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*array_function_type)(  ) const;
            array_function_type array_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::array );
            
            ChainPropertyProperty_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::assertCanConvert
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef void ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*assertCanConvert_function_type)( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::assertCanConvert );
            
            ChainPropertyProperty_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::assignFrom
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef void ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*assignFrom_function_type)( ::SireMol::ChainProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::assignFrom );
            
            ChainPropertyProperty_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::at
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > const & ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*at_function_type)( ::SireMol::ChainIdx const & ) const;
            at_function_type at_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::at );
            
            ChainPropertyProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("chainidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::at
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > const & ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::at );
            
            ChainPropertyProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::canConvert
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef bool ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*canConvert_function_type)( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::canConvert );
            
            ChainPropertyProperty_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::count
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef int ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::count );
            
            ChainPropertyProperty_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::fromVariant
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > ( *fromVariant_function_type )( ::SireMol::ChainProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::fromVariant );
            
            ChainPropertyProperty_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::get
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > const & ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*get_function_type)( ::SireMol::ChainIdx const & ) const;
            get_function_type get_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::get );
            
            ChainPropertyProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("chainidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::get
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > const & ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*get_function_type)( int ) const;
            get_function_type get_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::get );
            
            ChainPropertyProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::getAsProperty
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*getAsProperty_function_type)( ::SireMol::ChainIdx const & ) const;
            getAsProperty_function_type getAsProperty_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::getAsProperty );
            
            ChainPropertyProperty_exposer.def( 
                "getAsProperty"
                , getAsProperty_function_value
                , ( bp::arg("idx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::getAsVariant
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::QVariant ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*getAsVariant_function_type)( ::SireMol::ChainIdx const & ) const;
            getAsVariant_function_type getAsVariant_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::getAsVariant );
            
            ChainPropertyProperty_exposer.def( 
                "getAsVariant"
                , getAsVariant_function_value
                , ( bp::arg("idx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::isCompatibleWith
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef bool ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::isCompatibleWith );
            
            ChainPropertyProperty_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::isEmpty
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef bool ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::isEmpty );
            
            ChainPropertyProperty_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::nChains
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef int ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*nChains_function_type)(  ) const;
            nChains_function_type nChains_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::nChains );
            
            ChainPropertyProperty_exposer.def( 
                "nChains"
                , nChains_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ChainPropertyProperty_exposer.def( bp::self != bp::self );
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::operator=
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > & ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*assign_function_type)( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > const & ) ;
            assign_function_type assign_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::operator= );
            
            ChainPropertyProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ChainPropertyProperty_exposer.def( bp::self == bp::self );
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::operator[]
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > const & ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*__getitem___function_type)( ::SireMol::ChainIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::operator[] );
            
            ChainPropertyProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("chainidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::operator[]
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > const & ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::operator[] );
            
            ChainPropertyProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::operator[]
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::QList< SireBase::PropPtr< SireBase::Property > > ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*__getitem___function_type)( ::QList< long long > const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::operator[] );
            
            ChainPropertyProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idxs") )
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::operator[]
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::QList< SireBase::PropPtr< SireBase::Property > > ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*__getitem___function_type)( ::SireBase::Slice const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::operator[] );
            
            ChainPropertyProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::set
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > & ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*set_function_type)( ::SireMol::ChainIdx,::SireBase::PropPtr< SireBase::Property > const & ) ;
            set_function_type set_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::set );
            
            ChainPropertyProperty_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("chainidx"), bp::arg("value") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::size
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef int ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::size );
            
            ChainPropertyProperty_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::toString
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::QString ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::toString );
            
            ChainPropertyProperty_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::toVariant
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireMol::ChainProperty< QVariant > ( ::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::toVariant );
            
            ChainPropertyProperty_exposer.def( 
                "toVariant"
                , toVariant_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::typeName
        
            typedef SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ChainProperty< SireBase::PropPtr< SireBase::Property > >::typeName );
            
            ChainPropertyProperty_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ChainPropertyProperty_exposer.staticmethod( "fromVariant" );
        ChainPropertyProperty_exposer.staticmethod( "typeName" );
        ChainPropertyProperty_exposer.def( "__copy__", &__copy__);
        ChainPropertyProperty_exposer.def( "__deepcopy__", &__copy__);
        ChainPropertyProperty_exposer.def( "clone", &__copy__);
        ChainPropertyProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ChainProperty<SireBase::PropPtr<SireBase::Property> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainPropertyProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ChainProperty<SireBase::PropPtr<SireBase::Property> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainPropertyProperty_exposer.def_pickle(sire_pickle_suite< ::SireMol::ChainProperty<SireBase::PropPtr<SireBase::Property> > >());
        ChainPropertyProperty_exposer.def( "__str__", &__str__< ::SireMol::ChainProperty<SireBase::PropPtr<SireBase::Property> > > );
        ChainPropertyProperty_exposer.def( "__repr__", &__str__< ::SireMol::ChainProperty<SireBase::PropPtr<SireBase::Property> > > );
        ChainPropertyProperty_exposer.def( "__len__", &__len_size< ::SireMol::ChainProperty<SireBase::PropPtr<SireBase::Property> > > );
    }

}
