// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomPairs_LJScaleFactor_.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/parallel.h"

#include "SireMol/moleculeinfo.h"

#include "SireStream/datastream.h"

#include "cljnbpairs.h"

#include "cljnbpairs.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AtomPairs_LJScaleFactor__class(){

    { //::SireMM::AtomPairs< SireMM::LJScaleFactor >
        typedef bp::class_< SireMM::AtomPairs< SireMM::LJScaleFactor >, bp::bases< SireMol::MoleculeProperty, SireMol::MolViewProperty, SireBase::Property >, boost::noncopyable > AtomPairs_LJScaleFactor__exposer_t;
        AtomPairs_LJScaleFactor__exposer_t AtomPairs_LJScaleFactor__exposer = AtomPairs_LJScaleFactor__exposer_t( "AtomPairs_LJScaleFactor_", "", bp::no_init );
        bp::scope AtomPairs_LJScaleFactor__scope( AtomPairs_LJScaleFactor__exposer );
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::LJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*get_function_type)( ::SireMol::CGAtomIdx const & ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::get );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("atm0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::LJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*get_function_type)( ::SireMol::CGAtomIdx const &,::SireMol::CGAtomIdx const & ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::get );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("atm0"), bp::arg("atm1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::LJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*get_function_type)( ::SireMol::CGIdx ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::get );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgid0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::LJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*get_function_type)( ::SireMol::CGIdx,::SireMol::CGIdx ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::get );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgid0"), bp::arg("cgid1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::LJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*get_function_type)( ::SireMol::AtomID const & ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::get );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("atm0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::LJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*get_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::get );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("atm0"), bp::arg("atm1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::LJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*get_function_type)( ::SireMol::CGID const & ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::get );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgid0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::LJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*get_function_type)( ::SireMol::CGID const &,::SireMol::CGID const & ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::get );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgid0"), bp::arg("cgid1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::info
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMol::MoleculeInfoData const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*info_function_type)(  ) const;
            info_function_type info_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::info );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "info"
                , info_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::isCompatibleWith
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef bool ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::isCompatibleWith );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::isEmpty
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef bool ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::isEmpty );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::nAtoms
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef int ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::nAtoms );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::nGroups
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef int ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*nGroups_function_type)(  ) const;
            nGroups_function_type nGroups_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::nGroups );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "nGroups"
                , nGroups_function_value
                , "" );
        
        }
        AtomPairs_LJScaleFactor__exposer.def( bp::self != bp::self );
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::LJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*__call___function_type)( ::SireMol::CGAtomIdx const & ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator() );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("atm0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::LJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*__call___function_type)( ::SireMol::CGAtomIdx const &,::SireMol::CGAtomIdx const & ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator() );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("atm0"), bp::arg("atm1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::LJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*__call___function_type)( ::SireMol::CGIdx ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator() );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("cgid0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::LJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*__call___function_type)( ::SireMol::CGIdx,::SireMol::CGIdx ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator() );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("cgid0"), bp::arg("cgid1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::LJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*__call___function_type)( ::SireMol::AtomID const & ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator() );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("atm0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::LJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*__call___function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator() );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("atm0"), bp::arg("atm1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::LJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*__call___function_type)( ::SireMol::CGID const & ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator() );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("cgid0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::LJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*__call___function_type)( ::SireMol::CGID const &,::SireMol::CGID const & ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator() );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("cgid0"), bp::arg("cgid1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator=
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef ::SireMM::AtomPairs< SireMM::LJScaleFactor > & ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*assign_function_type)( ::SireMM::AtomPairs< SireMM::LJScaleFactor > const & ) ;
            assign_function_type assign_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::operator= );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AtomPairs_LJScaleFactor__exposer.def( bp::self == bp::self );
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::reserve
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*reserve_function_type)( int,int ) ;
            reserve_function_type reserve_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::reserve );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "reserve"
                , reserve_function_value
                , ( bp::arg("dim_x"), bp::arg("dim_y") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::set
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*set_function_type)( ::SireMol::CGAtomIdx const &,::SireMM::LJScaleFactor const & ) ;
            set_function_type set_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::set );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("atm0"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::set
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*set_function_type)( ::SireMol::CGAtomIdx const &,::SireMol::CGAtomIdx const &,::SireMM::LJScaleFactor const & ) ;
            set_function_type set_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::set );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("atm0"), bp::arg("atm1"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::set
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*set_function_type)( ::SireMol::AtomID const &,::SireMM::LJScaleFactor const & ) ;
            set_function_type set_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::set );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("atm0"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::set
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*set_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMM::LJScaleFactor const & ) ;
            set_function_type set_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::set );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("atm0"), bp::arg("atm1"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::setAll
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*setAll_function_type)( ::SireMM::LJScaleFactor const & ) ;
            setAll_function_type setAll_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::setAll );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::setAll
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*setAll_function_type)( ::SireMol::CGIdx,::SireMM::LJScaleFactor const & ) ;
            setAll_function_type setAll_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::setAll );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("cgid0"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::setAll
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*setAll_function_type)( ::SireMol::CGIdx,::SireMol::CGIdx,::SireMM::LJScaleFactor const & ) ;
            setAll_function_type setAll_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::setAll );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("cgid0"), bp::arg("cgid1"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::setAll
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*setAll_function_type)( ::SireMol::CGID const &,::SireMM::LJScaleFactor const & ) ;
            setAll_function_type setAll_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::setAll );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("cgid0"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::setAll
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*setAll_function_type)( ::SireMol::CGID const &,::SireMol::CGID const &,::SireMM::LJScaleFactor const & ) ;
            setAll_function_type setAll_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::setAll );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("cgid0"), bp::arg("cgid1"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::LJScaleFactor >::squeeze
        
            typedef SireMM::AtomPairs< SireMM::LJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::LJScaleFactor >::*squeeze_function_type)(  ) ;
            squeeze_function_type squeeze_function_value( &::SireMM::AtomPairs< SireMM::LJScaleFactor >::squeeze );
            
            AtomPairs_LJScaleFactor__exposer.def( 
                "squeeze"
                , squeeze_function_value
                , "" );
        
        }
        AtomPairs_LJScaleFactor__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::AtomPairs<SireMM::LJScaleFactor> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomPairs_LJScaleFactor__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::AtomPairs<SireMM::LJScaleFactor> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomPairs_LJScaleFactor__exposer.def( "__setstate__", &__setstate__base64< ::SireMM::AtomPairs<SireMM::LJScaleFactor> > );
        AtomPairs_LJScaleFactor__exposer.def( "__getstate__", &__getstate__base64< ::SireMM::AtomPairs<SireMM::LJScaleFactor> > );
        AtomPairs_LJScaleFactor__exposer.def( "__str__", &__str__< ::SireMM::AtomPairs<SireMM::LJScaleFactor> > );
        AtomPairs_LJScaleFactor__exposer.def( "__repr__", &__str__< ::SireMM::AtomPairs<SireMM::LJScaleFactor> > );
    }

}
