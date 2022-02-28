// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomPairs_CLJScaleFactor_.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/parallel.h"

#include "SireMol/moleculeinfo.h"

#include "SireStream/datastream.h"

#include "cljnbpairs.h"

#include "cljnbpairs.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AtomPairs_CLJScaleFactor__class(){

    { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >
        typedef bp::class_< SireMM::AtomPairs< SireMM::CLJScaleFactor >, bp::bases< SireMol::MoleculeProperty, SireMol::MolViewProperty, SireBase::Property >, boost::noncopyable > AtomPairs_CLJScaleFactor__exposer_t;
        AtomPairs_CLJScaleFactor__exposer_t AtomPairs_CLJScaleFactor__exposer = AtomPairs_CLJScaleFactor__exposer_t( "AtomPairs_CLJScaleFactor_", "", bp::no_init );
        bp::scope AtomPairs_CLJScaleFactor__scope( AtomPairs_CLJScaleFactor__exposer );
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CLJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*get_function_type)( ::SireMol::CGAtomIdx const & ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("atm0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CLJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*get_function_type)( ::SireMol::CGAtomIdx const &,::SireMol::CGAtomIdx const & ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("atm0"), bp::arg("atm1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::CLJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*get_function_type)( ::SireMol::CGIdx ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgid0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::CLJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*get_function_type)( ::SireMol::CGIdx,::SireMol::CGIdx ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgid0"), bp::arg("cgid1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CLJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*get_function_type)( ::SireMol::AtomID const & ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("atm0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CLJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*get_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("atm0"), bp::arg("atm1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::CLJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*get_function_type)( ::SireMol::CGID const & ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgid0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::CLJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*get_function_type)( ::SireMol::CGID const &,::SireMol::CGID const & ) const;
            get_function_type get_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::get );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgid0"), bp::arg("cgid1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::info
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMol::MoleculeInfoData const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*info_function_type)(  ) const;
            info_function_type info_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::info );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "info"
                , info_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::isCompatibleWith
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef bool ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::isCompatibleWith );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::isEmpty
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef bool ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::isEmpty );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::nAtoms
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef int ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::nAtoms );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::nGroups
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef int ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*nGroups_function_type)(  ) const;
            nGroups_function_type nGroups_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::nGroups );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "nGroups"
                , nGroups_function_value
                , "" );
        
        }
        AtomPairs_CLJScaleFactor__exposer.def( bp::self != bp::self );
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CLJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*__call___function_type)( ::SireMol::CGAtomIdx const & ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator() );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("atm0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CLJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*__call___function_type)( ::SireMol::CGAtomIdx const &,::SireMol::CGAtomIdx const & ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator() );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("atm0"), bp::arg("atm1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::CLJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*__call___function_type)( ::SireMol::CGIdx ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator() );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("cgid0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::CLJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*__call___function_type)( ::SireMol::CGIdx,::SireMol::CGIdx ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator() );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("cgid0"), bp::arg("cgid1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CLJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*__call___function_type)( ::SireMol::AtomID const & ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator() );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("atm0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CLJScaleFactor const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*__call___function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator() );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("atm0"), bp::arg("atm1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::CLJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*__call___function_type)( ::SireMol::CGID const & ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator() );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("cgid0") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator()
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::CGAtomPairs< SireMM::CLJScaleFactor > const & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*__call___function_type)( ::SireMol::CGID const &,::SireMol::CGID const & ) const;
            __call___function_type __call___function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator() );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("cgid0"), bp::arg("cgid1") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator=
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef ::SireMM::AtomPairs< SireMM::CLJScaleFactor > & ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*assign_function_type)( ::SireMM::AtomPairs< SireMM::CLJScaleFactor > const & ) ;
            assign_function_type assign_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::operator= );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AtomPairs_CLJScaleFactor__exposer.def( bp::self == bp::self );
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::reserve
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*reserve_function_type)( int,int ) ;
            reserve_function_type reserve_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::reserve );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "reserve"
                , reserve_function_value
                , ( bp::arg("dim_x"), bp::arg("dim_y") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::set
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*set_function_type)( ::SireMol::CGAtomIdx const &,::SireMM::CLJScaleFactor const & ) ;
            set_function_type set_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::set );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("atm0"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::set
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*set_function_type)( ::SireMol::CGAtomIdx const &,::SireMol::CGAtomIdx const &,::SireMM::CLJScaleFactor const & ) ;
            set_function_type set_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::set );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("atm0"), bp::arg("atm1"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::set
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*set_function_type)( ::SireMol::AtomID const &,::SireMM::CLJScaleFactor const & ) ;
            set_function_type set_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::set );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("atm0"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::set
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*set_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMM::CLJScaleFactor const & ) ;
            set_function_type set_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::set );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("atm0"), bp::arg("atm1"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::setAll
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*setAll_function_type)( ::SireMM::CLJScaleFactor const & ) ;
            setAll_function_type setAll_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::setAll );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::setAll
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*setAll_function_type)( ::SireMol::CGIdx,::SireMM::CLJScaleFactor const & ) ;
            setAll_function_type setAll_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::setAll );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("cgid0"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::setAll
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*setAll_function_type)( ::SireMol::CGIdx,::SireMol::CGIdx,::SireMM::CLJScaleFactor const & ) ;
            setAll_function_type setAll_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::setAll );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("cgid0"), bp::arg("cgid1"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::setAll
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*setAll_function_type)( ::SireMol::CGID const &,::SireMM::CLJScaleFactor const & ) ;
            setAll_function_type setAll_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::setAll );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("cgid0"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::setAll
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*setAll_function_type)( ::SireMol::CGID const &,::SireMol::CGID const &,::SireMM::CLJScaleFactor const & ) ;
            setAll_function_type setAll_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::setAll );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("cgid0"), bp::arg("cgid1"), bp::arg("value") )
                , "" );
        
        }
        { //::SireMM::AtomPairs< SireMM::CLJScaleFactor >::squeeze
        
            typedef SireMM::AtomPairs< SireMM::CLJScaleFactor > exported_class_t;
            typedef void ( ::SireMM::AtomPairs< SireMM::CLJScaleFactor >::*squeeze_function_type)(  ) ;
            squeeze_function_type squeeze_function_value( &::SireMM::AtomPairs< SireMM::CLJScaleFactor >::squeeze );
            
            AtomPairs_CLJScaleFactor__exposer.def( 
                "squeeze"
                , squeeze_function_value
                , "" );
        
        }
        AtomPairs_CLJScaleFactor__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::AtomPairs<SireMM::CLJScaleFactor> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomPairs_CLJScaleFactor__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::AtomPairs<SireMM::CLJScaleFactor> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomPairs_CLJScaleFactor__exposer.def_pickle(sire_pickle_suite< ::SireMM::AtomPairs<SireMM::CLJScaleFactor> >());
        AtomPairs_CLJScaleFactor__exposer.def( "__str__", &__str__< ::SireMM::AtomPairs<SireMM::CLJScaleFactor> > );
        AtomPairs_CLJScaleFactor__exposer.def( "__repr__", &__str__< ::SireMM::AtomPairs<SireMM::CLJScaleFactor> > );
    }

}
