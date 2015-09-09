// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "MultiDouble.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "multidouble.h"

#include <QDebug>

#include <QStringList>

#include "multidouble.h"

#include "multifloat.h"

#include "multiint.h"

SireMaths::MultiDouble __copy__(const SireMaths::MultiDouble &other){ return SireMaths::MultiDouble(other); }

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_MultiDouble_class(){

    { //::SireMaths::MultiDouble
        typedef bp::class_< SireMaths::MultiDouble > MultiDouble_exposer_t;
        MultiDouble_exposer_t MultiDouble_exposer = MultiDouble_exposer_t( "MultiDouble", bp::init< >() );
        bp::scope MultiDouble_scope( MultiDouble_exposer );
        MultiDouble_exposer.def( bp::init< double >(( bp::arg("value") )) );
        MultiDouble_exposer.def( bp::init< double const *, int >(( bp::arg("array"), bp::arg("size") )) );
        MultiDouble_exposer.def( bp::init< QVector< float > const & >(( bp::arg("array") )) );
        MultiDouble_exposer.def( bp::init< QVector< double > const & >(( bp::arg("array") )) );
        MultiDouble_exposer.def( bp::init< SireMaths::MultiFloat const & >(( bp::arg("other") )) );
        MultiDouble_exposer.def( bp::init< SireMaths::MultiDouble const & >(( bp::arg("other") )) );
        { //::SireMaths::MultiDouble::at
        
            typedef double ( ::SireMaths::MultiDouble::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMaths::MultiDouble::at );
            
            MultiDouble_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMaths::MultiDouble::compareEqual
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*compareEqual_function_type)( ::SireMaths::MultiDouble const & ) const;
            compareEqual_function_type compareEqual_function_value( &::SireMaths::MultiDouble::compareEqual );
            
            MultiDouble_exposer.def( 
                "compareEqual"
                , compareEqual_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::MultiDouble::compareGreater
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*compareGreater_function_type)( ::SireMaths::MultiDouble const & ) const;
            compareGreater_function_type compareGreater_function_value( &::SireMaths::MultiDouble::compareGreater );
            
            MultiDouble_exposer.def( 
                "compareGreater"
                , compareGreater_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::MultiDouble::compareGreaterEqual
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*compareGreaterEqual_function_type)( ::SireMaths::MultiDouble const & ) const;
            compareGreaterEqual_function_type compareGreaterEqual_function_value( &::SireMaths::MultiDouble::compareGreaterEqual );
            
            MultiDouble_exposer.def( 
                "compareGreaterEqual"
                , compareGreaterEqual_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::MultiDouble::compareLess
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*compareLess_function_type)( ::SireMaths::MultiDouble const & ) const;
            compareLess_function_type compareLess_function_value( &::SireMaths::MultiDouble::compareLess );
            
            MultiDouble_exposer.def( 
                "compareLess"
                , compareLess_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::MultiDouble::compareLessEqual
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*compareLessEqual_function_type)( ::SireMaths::MultiDouble const & ) const;
            compareLessEqual_function_type compareLessEqual_function_value( &::SireMaths::MultiDouble::compareLessEqual );
            
            MultiDouble_exposer.def( 
                "compareLessEqual"
                , compareLessEqual_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::MultiDouble::compareNotEqual
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*compareNotEqual_function_type)( ::SireMaths::MultiDouble const & ) const;
            compareNotEqual_function_type compareNotEqual_function_value( &::SireMaths::MultiDouble::compareNotEqual );
            
            MultiDouble_exposer.def( 
                "compareNotEqual"
                , compareNotEqual_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::MultiDouble::count
        
            typedef int ( *count_function_type )(  );
            count_function_type count_function_value( &::SireMaths::MultiDouble::count );
            
            MultiDouble_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireMaths::MultiDouble::doubleSum
        
            typedef double ( ::SireMaths::MultiDouble::*doubleSum_function_type)(  ) const;
            doubleSum_function_type doubleSum_function_value( &::SireMaths::MultiDouble::doubleSum );
            
            MultiDouble_exposer.def( 
                "doubleSum"
                , doubleSum_function_value );
        
        }
        { //::SireMaths::MultiDouble::fromArray
        
            typedef ::QVector< SireMaths::MultiDouble > ( *fromArray_function_type )( ::QVector< double > const & );
            fromArray_function_type fromArray_function_value( &::SireMaths::MultiDouble::fromArray );
            
            MultiDouble_exposer.def( 
                "fromArray"
                , fromArray_function_value
                , ( bp::arg("array") ) );
        
        }
        { //::SireMaths::MultiDouble::fromArray
        
            typedef ::QVector< SireMaths::MultiDouble > ( *fromArray_function_type )( ::QVector< float > const & );
            fromArray_function_type fromArray_function_value( &::SireMaths::MultiDouble::fromArray );
            
            MultiDouble_exposer.def( 
                "fromArray"
                , fromArray_function_value
                , ( bp::arg("array") ) );
        
        }
        { //::SireMaths::MultiDouble::fromArray
        
            typedef ::QVector< SireMaths::MultiDouble > ( *fromArray_function_type )( double const *,int );
            fromArray_function_type fromArray_function_value( &::SireMaths::MultiDouble::fromArray );
            
            MultiDouble_exposer.def( 
                "fromArray"
                , fromArray_function_value
                , ( bp::arg("array"), bp::arg("size") ) );
        
        }
        { //::SireMaths::MultiDouble::fromArray
        
            typedef ::QVector< SireMaths::MultiDouble > ( *fromArray_function_type )( float const *,int );
            fromArray_function_type fromArray_function_value( &::SireMaths::MultiDouble::fromArray );
            
            MultiDouble_exposer.def( 
                "fromArray"
                , fromArray_function_value
                , ( bp::arg("array"), bp::arg("size") ) );
        
        }
        { //::SireMaths::MultiDouble::get
        
            typedef double ( ::SireMaths::MultiDouble::*get_function_type)( int ) const;
            get_function_type get_function_value( &::SireMaths::MultiDouble::get );
            
            MultiDouble_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMaths::MultiDouble::getitem
        
            typedef double ( ::SireMaths::MultiDouble::*getitem_function_type)( int ) const;
            getitem_function_type getitem_function_value( &::SireMaths::MultiDouble::getitem );
            
            MultiDouble_exposer.def( 
                "getitem"
                , getitem_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMaths::MultiDouble::isAligned
        
            typedef bool ( ::SireMaths::MultiDouble::*isAligned_function_type)(  ) const;
            isAligned_function_type isAligned_function_value( &::SireMaths::MultiDouble::isAligned );
            
            MultiDouble_exposer.def( 
                "isAligned"
                , isAligned_function_value );
        
        }
        { //::SireMaths::MultiDouble::logicalAnd
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*logicalAnd_function_type)( ::SireMaths::MultiDouble const & ) const;
            logicalAnd_function_type logicalAnd_function_value( &::SireMaths::MultiDouble::logicalAnd );
            
            MultiDouble_exposer.def( 
                "logicalAnd"
                , logicalAnd_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::MultiDouble::logicalAndNot
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*logicalAndNot_function_type)( ::SireMaths::MultiDouble const & ) const;
            logicalAndNot_function_type logicalAndNot_function_value( &::SireMaths::MultiDouble::logicalAndNot );
            
            MultiDouble_exposer.def( 
                "logicalAndNot"
                , logicalAndNot_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::MultiDouble::logicalNot
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*logicalNot_function_type)(  ) const;
            logicalNot_function_type logicalNot_function_value( &::SireMaths::MultiDouble::logicalNot );
            
            MultiDouble_exposer.def( 
                "logicalNot"
                , logicalNot_function_value );
        
        }
        { //::SireMaths::MultiDouble::logicalOr
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*logicalOr_function_type)( ::SireMaths::MultiDouble const & ) const;
            logicalOr_function_type logicalOr_function_value( &::SireMaths::MultiDouble::logicalOr );
            
            MultiDouble_exposer.def( 
                "logicalOr"
                , logicalOr_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::MultiDouble::logicalXor
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*logicalXor_function_type)( ::SireMaths::MultiDouble const & ) const;
            logicalXor_function_type logicalXor_function_value( &::SireMaths::MultiDouble::logicalXor );
            
            MultiDouble_exposer.def( 
                "logicalXor"
                , logicalXor_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::MultiDouble::max
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*max_function_type)( ::SireMaths::MultiDouble const & ) const;
            max_function_type max_function_value( &::SireMaths::MultiDouble::max );
            
            MultiDouble_exposer.def( 
                "max"
                , max_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::MultiDouble::min
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*min_function_type)( ::SireMaths::MultiDouble const & ) const;
            min_function_type min_function_value( &::SireMaths::MultiDouble::min );
            
            MultiDouble_exposer.def( 
                "min"
                , min_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::MultiDouble::multiplyAdd
        
            typedef ::SireMaths::MultiDouble & ( ::SireMaths::MultiDouble::*multiplyAdd_function_type)( ::SireMaths::MultiDouble const &,::SireMaths::MultiDouble const & ) ;
            multiplyAdd_function_type multiplyAdd_function_value( &::SireMaths::MultiDouble::multiplyAdd );
            
            MultiDouble_exposer.def( 
                "multiplyAdd"
                , multiplyAdd_function_value
                , ( bp::arg("val0"), bp::arg("val1") )
                , bp::return_self< >() );
        
        }
        MultiDouble_exposer.def( !bp::self );
        MultiDouble_exposer.def( bp::self != bp::self );
        MultiDouble_exposer.def( bp::self & bp::self );
        MultiDouble_exposer.def( bp::self * bp::self );
        MultiDouble_exposer.def( bp::self + bp::self );
        MultiDouble_exposer.def( -bp::self );
        MultiDouble_exposer.def( bp::self - bp::self );
        MultiDouble_exposer.def( bp::self / bp::self );
        MultiDouble_exposer.def( bp::self < bp::self );
        MultiDouble_exposer.def( bp::self <= bp::self );
        { //::SireMaths::MultiDouble::operator=
        
            typedef ::SireMaths::MultiDouble & ( ::SireMaths::MultiDouble::*assign_function_type)( ::SireMaths::MultiDouble const & ) ;
            assign_function_type assign_function_value( &::SireMaths::MultiDouble::operator= );
            
            MultiDouble_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireMaths::MultiDouble::operator=
        
            typedef ::SireMaths::MultiDouble & ( ::SireMaths::MultiDouble::*assign_function_type)( ::SireMaths::MultiFloat const & ) ;
            assign_function_type assign_function_value( &::SireMaths::MultiDouble::operator= );
            
            MultiDouble_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireMaths::MultiDouble::operator=
        
            typedef ::SireMaths::MultiDouble & ( ::SireMaths::MultiDouble::*assign_function_type)( double ) ;
            assign_function_type assign_function_value( &::SireMaths::MultiDouble::operator= );
            
            MultiDouble_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("value") )
                , bp::return_self< >() );
        
        }
        MultiDouble_exposer.def( bp::self == bp::self );
        MultiDouble_exposer.def( bp::self > bp::self );
        MultiDouble_exposer.def( bp::self >= bp::self );
        { //::SireMaths::MultiDouble::operator[]
        
            typedef double ( ::SireMaths::MultiDouble::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMaths::MultiDouble::operator[] );
            
            MultiDouble_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") ) );
        
        }
        MultiDouble_exposer.def( bp::self ^ bp::self );
        MultiDouble_exposer.def( bp::self | bp::self );
        { //::SireMaths::MultiDouble::reciprocal
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*reciprocal_function_type)(  ) const;
            reciprocal_function_type reciprocal_function_value( &::SireMaths::MultiDouble::reciprocal );
            
            MultiDouble_exposer.def( 
                "reciprocal"
                , reciprocal_function_value );
        
        }
        { //::SireMaths::MultiDouble::rotate
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*rotate_function_type)(  ) const;
            rotate_function_type rotate_function_value( &::SireMaths::MultiDouble::rotate );
            
            MultiDouble_exposer.def( 
                "rotate"
                , rotate_function_value );
        
        }
        { //::SireMaths::MultiDouble::rsqrt
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*rsqrt_function_type)(  ) const;
            rsqrt_function_type rsqrt_function_value( &::SireMaths::MultiDouble::rsqrt );
            
            MultiDouble_exposer.def( 
                "rsqrt"
                , rsqrt_function_value );
        
        }
        { //::SireMaths::MultiDouble::set
        
            typedef void ( ::SireMaths::MultiDouble::*set_function_type)( int,double ) ;
            set_function_type set_function_value( &::SireMaths::MultiDouble::set );
            
            MultiDouble_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("value") ) );
        
        }
        { //::SireMaths::MultiDouble::size
        
            typedef int ( *size_function_type )(  );
            size_function_type size_function_value( &::SireMaths::MultiDouble::size );
            
            MultiDouble_exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::SireMaths::MultiDouble::sqrt
        
            typedef ::SireMaths::MultiDouble ( ::SireMaths::MultiDouble::*sqrt_function_type)(  ) const;
            sqrt_function_type sqrt_function_value( &::SireMaths::MultiDouble::sqrt );
            
            MultiDouble_exposer.def( 
                "sqrt"
                , sqrt_function_value );
        
        }
        { //::SireMaths::MultiDouble::sum
        
            typedef double ( ::SireMaths::MultiDouble::*sum_function_type)(  ) const;
            sum_function_type sum_function_value( &::SireMaths::MultiDouble::sum );
            
            MultiDouble_exposer.def( 
                "sum"
                , sum_function_value );
        
        }
        { //::SireMaths::MultiDouble::toArray
        
            typedef ::QVector< double > ( *toArray_function_type )( ::QVector< SireMaths::MultiDouble > const & );
            toArray_function_type toArray_function_value( &::SireMaths::MultiDouble::toArray );
            
            MultiDouble_exposer.def( 
                "toArray"
                , toArray_function_value
                , ( bp::arg("array") ) );
        
        }
        { //::SireMaths::MultiDouble::toBinaryString
        
            typedef ::QString ( ::SireMaths::MultiDouble::*toBinaryString_function_type)(  ) const;
            toBinaryString_function_type toBinaryString_function_value( &::SireMaths::MultiDouble::toBinaryString );
            
            MultiDouble_exposer.def( 
                "toBinaryString"
                , toBinaryString_function_value );
        
        }
        { //::SireMaths::MultiDouble::toDoubleArray
        
            typedef ::QVector< double > ( *toDoubleArray_function_type )( ::QVector< SireMaths::MultiDouble > const & );
            toDoubleArray_function_type toDoubleArray_function_value( &::SireMaths::MultiDouble::toDoubleArray );
            
            MultiDouble_exposer.def( 
                "toDoubleArray"
                , toDoubleArray_function_value
                , ( bp::arg("array") ) );
        
        }
        { //::SireMaths::MultiDouble::toString
        
            typedef ::QString ( ::SireMaths::MultiDouble::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMaths::MultiDouble::toString );
            
            MultiDouble_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMaths::MultiDouble::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::MultiDouble::typeName );
            
            MultiDouble_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMaths::MultiDouble::what
        
            typedef char const * ( ::SireMaths::MultiDouble::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMaths::MultiDouble::what );
            
            MultiDouble_exposer.def( 
                "what"
                , what_function_value );
        
        }
        MultiDouble_exposer.staticmethod( "count" );
        MultiDouble_exposer.staticmethod( "fromArray" );
        MultiDouble_exposer.staticmethod( "size" );
        MultiDouble_exposer.staticmethod( "toArray" );
        MultiDouble_exposer.staticmethod( "toDoubleArray" );
        MultiDouble_exposer.staticmethod( "typeName" );
        MultiDouble_exposer.def( "__copy__", &__copy__);
        MultiDouble_exposer.def( "__deepcopy__", &__copy__);
        MultiDouble_exposer.def( "clone", &__copy__);
        MultiDouble_exposer.def( "__str__", &__str__< ::SireMaths::MultiDouble > );
        MultiDouble_exposer.def( "__repr__", &__str__< ::SireMaths::MultiDouble > );
        MultiDouble_exposer.def( "__len__", &__len_size< ::SireMaths::MultiDouble > );
        MultiDouble_exposer.def( "__getitem__", &::SireMaths::MultiDouble::getitem );
    }

}
