// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ZMatrixCoordsLine.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMol/angleid.h"

#include "SireMol/atommatcher.h"

#include "SireMol/bondid.h"

#include "SireMol/dihedralid.h"

#include "SireMol/errors.h"

#include "SireMol/molecule.h"

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "SireMove/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/convert.h"

#include "SireUnits/units.h"

#include "zmatrix.h"

#include <QDebug>

#include <QElapsedTimer>

#include <QTime>

#include "zmatrix.h"

SireMove::ZMatrixCoordsLine __copy__(const SireMove::ZMatrixCoordsLine &other){ return SireMove::ZMatrixCoordsLine(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ZMatrixCoordsLine_class(){

    { //::SireMove::ZMatrixCoordsLine
        typedef bp::class_< SireMove::ZMatrixCoordsLine, bp::bases< SireMove::ZMatrixLine > > ZMatrixCoordsLine_exposer_t;
        ZMatrixCoordsLine_exposer_t ZMatrixCoordsLine_exposer = ZMatrixCoordsLine_exposer_t( "ZMatrixCoordsLine", "This class holds a z-matrix line that includes the\nsizes of the internal coordinates\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope ZMatrixCoordsLine_scope( ZMatrixCoordsLine_exposer );
        ZMatrixCoordsLine_exposer.def( bp::init< SireMove::ZMatrixLine const & >(( bp::arg("line") ), "Construct from the passed line (but with zero coordinates)") );
        ZMatrixCoordsLine_exposer.def( bp::init< SireMove::ZMatrixLine const &, SireUnits::Dimension::Length const &, SireUnits::Dimension::Angle const &, SireUnits::Dimension::Angle const & >(( bp::arg("line"), bp::arg("bond"), bp::arg("angle"), bp::arg("dihedral") ), "Construct from the passed line and coordinates") );
        ZMatrixCoordsLine_exposer.def( bp::init< SireMove::ZMatrixCoordsLine const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::ZMatrixCoordsLine::angleSize
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMove::ZMatrixCoordsLine::*angleSize_function_type)(  ) const;
            angleSize_function_type angleSize_function_value( &::SireMove::ZMatrixCoordsLine::angleSize );
            
            ZMatrixCoordsLine_exposer.def( 
                "angleSize"
                , angleSize_function_value
                , "Return the size of the angle" );
        
        }
        { //::SireMove::ZMatrixCoordsLine::bondLength
        
            typedef ::SireUnits::Dimension::Length ( ::SireMove::ZMatrixCoordsLine::*bondLength_function_type)(  ) const;
            bondLength_function_type bondLength_function_value( &::SireMove::ZMatrixCoordsLine::bondLength );
            
            ZMatrixCoordsLine_exposer.def( 
                "bondLength"
                , bondLength_function_value
                , "Return the length of the bond" );
        
        }
        { //::SireMove::ZMatrixCoordsLine::dihedralSize
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMove::ZMatrixCoordsLine::*dihedralSize_function_type)(  ) const;
            dihedralSize_function_type dihedralSize_function_value( &::SireMove::ZMatrixCoordsLine::dihedralSize );
            
            ZMatrixCoordsLine_exposer.def( 
                "dihedralSize"
                , dihedralSize_function_value
                , "Return the size of the dihedral" );
        
        }
        ZMatrixCoordsLine_exposer.def( bp::self != bp::self );
        { //::SireMove::ZMatrixCoordsLine::operator=
        
            typedef ::SireMove::ZMatrixCoordsLine & ( ::SireMove::ZMatrixCoordsLine::*assign_function_type)( ::SireMove::ZMatrixCoordsLine const & ) ;
            assign_function_type assign_function_value( &::SireMove::ZMatrixCoordsLine::operator= );
            
            ZMatrixCoordsLine_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ZMatrixCoordsLine_exposer.def( bp::self == bp::self );
        { //::SireMove::ZMatrixCoordsLine::setAngle
        
            typedef void ( ::SireMove::ZMatrixCoordsLine::*setAngle_function_type)( ::SireUnits::Dimension::Angle const & ) ;
            setAngle_function_type setAngle_function_value( &::SireMove::ZMatrixCoordsLine::setAngle );
            
            ZMatrixCoordsLine_exposer.def( 
                "setAngle"
                , setAngle_function_value
                , ( bp::arg("size") )
                , "Set the size of the angle" );
        
        }
        { //::SireMove::ZMatrixCoordsLine::setBond
        
            typedef void ( ::SireMove::ZMatrixCoordsLine::*setBond_function_type)( ::SireUnits::Dimension::Length const & ) ;
            setBond_function_type setBond_function_value( &::SireMove::ZMatrixCoordsLine::setBond );
            
            ZMatrixCoordsLine_exposer.def( 
                "setBond"
                , setBond_function_value
                , ( bp::arg("length") )
                , "Set the length of the bond" );
        
        }
        { //::SireMove::ZMatrixCoordsLine::setDihedral
        
            typedef void ( ::SireMove::ZMatrixCoordsLine::*setDihedral_function_type)( ::SireUnits::Dimension::Angle const & ) ;
            setDihedral_function_type setDihedral_function_value( &::SireMove::ZMatrixCoordsLine::setDihedral );
            
            ZMatrixCoordsLine_exposer.def( 
                "setDihedral"
                , setDihedral_function_value
                , ( bp::arg("size") )
                , "Set the size of the dihedral" );
        
        }
        { //::SireMove::ZMatrixCoordsLine::toString
        
            typedef ::QString ( ::SireMove::ZMatrixCoordsLine::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMove::ZMatrixCoordsLine::toString );
            
            ZMatrixCoordsLine_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation" );
        
        }
        { //::SireMove::ZMatrixCoordsLine::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::ZMatrixCoordsLine::typeName );
            
            ZMatrixCoordsLine_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMove::ZMatrixCoordsLine::what
        
            typedef char const * ( ::SireMove::ZMatrixCoordsLine::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMove::ZMatrixCoordsLine::what );
            
            ZMatrixCoordsLine_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        ZMatrixCoordsLine_exposer.staticmethod( "typeName" );
        ZMatrixCoordsLine_exposer.def( "__copy__", &__copy__);
        ZMatrixCoordsLine_exposer.def( "__deepcopy__", &__copy__);
        ZMatrixCoordsLine_exposer.def( "clone", &__copy__);
        ZMatrixCoordsLine_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::ZMatrixCoordsLine >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ZMatrixCoordsLine_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::ZMatrixCoordsLine >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ZMatrixCoordsLine_exposer.def( "__getstate_manages_dict__", true);
        ZMatrixCoordsLine_exposer.def( "__safe_for_unpickling__", true);
        ZMatrixCoordsLine_exposer.def( "__setstate__", &__setstate__base64< ::SireMove::ZMatrixCoordsLine > );
        ZMatrixCoordsLine_exposer.def( "__getstate__", &__getstate__base64< ::SireMove::ZMatrixCoordsLine > );
        ZMatrixCoordsLine_exposer.def( "__str__", &__str__< ::SireMove::ZMatrixCoordsLine > );
        ZMatrixCoordsLine_exposer.def( "__repr__", &__str__< ::SireMove::ZMatrixCoordsLine > );
    }

}
