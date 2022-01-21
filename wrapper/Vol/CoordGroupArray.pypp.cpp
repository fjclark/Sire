// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CoordGroupArray.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/quickcopy.hpp"

#include "SireError/errors.h"

#include "SireMaths/align.h"

#include "SireMaths/axisset.h"

#include "SireMaths/matrix.h"

#include "SireMaths/quaternion.h"

#include "SireMaths/rotate.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "coordgroup.h"

#include <QDebug>

#include "coordgroup.h"

SireVol::CoordGroupArray __copy__(const SireVol::CoordGroupArray &other){ return SireVol::CoordGroupArray(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_CoordGroupArray_class(){

    { //::SireVol::CoordGroupArray
        typedef bp::class_< SireVol::CoordGroupArray > CoordGroupArray_exposer_t;
        CoordGroupArray_exposer_t CoordGroupArray_exposer = CoordGroupArray_exposer_t( "CoordGroupArray", "This class holds an array of CoordGroups. While you could\nof course just use a QVector<CoordGroup>, this array\noptimises the memory layout of all of the CoordGroups\nso that they all lie contiguously along the same piece\nof memory (and indeed, all of the AABoxes are grouped\ntogether, while all of the coordinates are grouped together).\n\nThe memory packing means that this array is much more\nlimited than a QVector<CoordGroup>, i.e. you cant\nadd or remove CoordGroups from the array, and you cant\ndo anything to the contained CoordGroups except for\nchange their coordinates.\n\nThis class is really meant to be used as a fast container\nthat allow rapid iteration over all of the contained\nCoordGroups  coordinates\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope CoordGroupArray_scope( CoordGroupArray_exposer );
        CoordGroupArray_exposer.def( bp::init< SireVol::CoordGroup const & >(( bp::arg("cgroup") ), "Construct an array that holds just the passed CoordGroup") );
        CoordGroupArray_exposer.def( bp::init< QVector< QVector< SireMaths::Vector > > const & >(( bp::arg("points") ), "Construct from a double-vector") );
        CoordGroupArray_exposer.def( bp::init< QVector< SireVol::CoordGroup > const & >(( bp::arg("cgroups") ), "Construct from an array of CoordGroups") );
        CoordGroupArray_exposer.def( bp::init< SireVol::CoordGroupArray const &, SireVol::CoordGroupArray const & >(( bp::arg("array0"), bp::arg("array1") ), "Construct from a pair of CoordGroupArrays") );
        CoordGroupArray_exposer.def( bp::init< SireVol::CoordGroupArray const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireVol::CoordGroupArray::aaBox
        
            typedef ::SireVol::AABox ( ::SireVol::CoordGroupArray::*aaBox_function_type)(  ) const;
            aaBox_function_type aaBox_function_value( &::SireVol::CoordGroupArray::aaBox );
            
            CoordGroupArray_exposer.def( 
                "aaBox"
                , aaBox_function_value
                , "Return an AABox that complete encompasses all of the CoordGroups\nin this array" );
        
        }
        { //::SireVol::CoordGroupArray::append
        
            typedef void ( ::SireVol::CoordGroupArray::*append_function_type)( ::SireVol::CoordGroup const & ) ;
            append_function_type append_function_value( &::SireVol::CoordGroupArray::append );
            
            CoordGroupArray_exposer.def( 
                "append"
                , append_function_value
                , ( bp::arg("cgroup") )
                , "Append the passed CoordGroup onto the end of this array" );
        
        }
        { //::SireVol::CoordGroupArray::append
        
            typedef void ( ::SireVol::CoordGroupArray::*append_function_type)( ::SireVol::CoordGroupArray const & ) ;
            append_function_type append_function_value( &::SireVol::CoordGroupArray::append );
            
            CoordGroupArray_exposer.def( 
                "append"
                , append_function_value
                , ( bp::arg("cgroups") )
                , "Append the passed CoordGroups onto the end of this array" );
        
        }
        { //::SireVol::CoordGroupArray::assertValidCoordGroup
        
            typedef void ( ::SireVol::CoordGroupArray::*assertValidCoordGroup_function_type)( ::quint32 ) const;
            assertValidCoordGroup_function_type assertValidCoordGroup_function_value( &::SireVol::CoordGroupArray::assertValidCoordGroup );
            
            CoordGroupArray_exposer.def( 
                "assertValidCoordGroup"
                , assertValidCoordGroup_function_value
                , ( bp::arg("i") )
                , "Assert that the index i points to a valid CoordGroup\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireVol::CoordGroupArray::assertValidCoordinate
        
            typedef void ( ::SireVol::CoordGroupArray::*assertValidCoordinate_function_type)( ::quint32 ) const;
            assertValidCoordinate_function_type assertValidCoordinate_function_value( &::SireVol::CoordGroupArray::assertValidCoordinate );
            
            CoordGroupArray_exposer.def( 
                "assertValidCoordinate"
                , assertValidCoordinate_function_value
                , ( bp::arg("i") )
                , "Assert that the index i points a valid coordinate\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireVol::CoordGroupArray::assertValidIndex
        
            typedef void ( ::SireVol::CoordGroupArray::*assertValidIndex_function_type)( ::quint32 ) const;
            assertValidIndex_function_type assertValidIndex_function_value( &::SireVol::CoordGroupArray::assertValidIndex );
            
            CoordGroupArray_exposer.def( 
                "assertValidIndex"
                , assertValidIndex_function_value
                , ( bp::arg("i") )
                , "Assert that the index i points to a valid CoordGroup\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireVol::CoordGroupArray::at
        
            typedef ::SireVol::CoordGroup const & ( ::SireVol::CoordGroupArray::*at_function_type)( ::quint32 ) const;
            at_function_type at_function_value( &::SireVol::CoordGroupArray::at );
            
            CoordGroupArray_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return a reference to the ith CoordGroup\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireVol::CoordGroupArray::changeFrame
        
            typedef void ( ::SireVol::CoordGroupArray::*changeFrame_function_type)( ::SireMaths::AxisSet const &,::SireMaths::AxisSet const & ) ;
            changeFrame_function_type changeFrame_function_value( &::SireVol::CoordGroupArray::changeFrame );
            
            CoordGroupArray_exposer.def( 
                "changeFrame"
                , changeFrame_function_value
                , ( bp::arg("from_frame"), bp::arg("to_frame") )
                , "Change all of the coordinates in this array from the\ncoordinate frame from_frame to the coordinate frame to_frame" );
        
        }
        { //::SireVol::CoordGroupArray::changeFrame
        
            typedef void ( ::SireVol::CoordGroupArray::*changeFrame_function_type)( ::quint32,::SireMaths::AxisSet const &,::SireMaths::AxisSet const & ) ;
            changeFrame_function_type changeFrame_function_value( &::SireVol::CoordGroupArray::changeFrame );
            
            CoordGroupArray_exposer.def( 
                "changeFrame"
                , changeFrame_function_value
                , ( bp::arg("i"), bp::arg("from_frame"), bp::arg("to_frame") )
                , "Change all of the coordinates in the ith CoordGroup from\nthe coordinate frame from_frame to the coordinate frame\nto_frame\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireVol::CoordGroupArray::count
        
            typedef int ( ::SireVol::CoordGroupArray::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireVol::CoordGroupArray::count );
            
            CoordGroupArray_exposer.def( 
                "count"
                , count_function_value
                , "Return the number of CoordGroups in this array" );
        
        }
        { //::SireVol::CoordGroupArray::isEmpty
        
            typedef bool ( ::SireVol::CoordGroupArray::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireVol::CoordGroupArray::isEmpty );
            
            CoordGroupArray_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , "Return whether or not this array is empty" );
        
        }
        { //::SireVol::CoordGroupArray::mapInto
        
            typedef void ( ::SireVol::CoordGroupArray::*mapInto_function_type)( ::SireMaths::AxisSet const & ) ;
            mapInto_function_type mapInto_function_value( &::SireVol::CoordGroupArray::mapInto );
            
            CoordGroupArray_exposer.def( 
                "mapInto"
                , mapInto_function_value
                , ( bp::arg("axes") )
                , "Map all of the coordinates in this array into the coordinate\nframe represented by axes" );
        
        }
        { //::SireVol::CoordGroupArray::mapInto
        
            typedef void ( ::SireVol::CoordGroupArray::*mapInto_function_type)( ::quint32,::SireMaths::AxisSet const & ) ;
            mapInto_function_type mapInto_function_value( &::SireVol::CoordGroupArray::mapInto );
            
            CoordGroupArray_exposer.def( 
                "mapInto"
                , mapInto_function_value
                , ( bp::arg("i"), bp::arg("axes") )
                , "Map all of the coordinates of the CoordGroup at index i\ninto the coordinate frame represented by axes\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireVol::CoordGroupArray::merge
        
            typedef ::SireVol::CoordGroup ( ::SireVol::CoordGroupArray::*merge_function_type)(  ) const;
            merge_function_type merge_function_value( &::SireVol::CoordGroupArray::merge );
            
            CoordGroupArray_exposer.def( 
                "merge"
                , merge_function_value
                , "Merge this array of CoordGroups back into a single CoordGroup" );
        
        }
        { //::SireVol::CoordGroupArray::nCoordGroups
        
            typedef int ( ::SireVol::CoordGroupArray::*nCoordGroups_function_type)(  ) const;
            nCoordGroups_function_type nCoordGroups_function_value( &::SireVol::CoordGroupArray::nCoordGroups );
            
            CoordGroupArray_exposer.def( 
                "nCoordGroups"
                , nCoordGroups_function_value
                , "Return the number of CoordGroups in this array" );
        
        }
        { //::SireVol::CoordGroupArray::nCoords
        
            typedef int ( ::SireVol::CoordGroupArray::*nCoords_function_type)(  ) const;
            nCoords_function_type nCoords_function_value( &::SireVol::CoordGroupArray::nCoords );
            
            CoordGroupArray_exposer.def( 
                "nCoords"
                , nCoords_function_value
                , "Return the number of coordinates in this array" );
        
        }
        CoordGroupArray_exposer.def( bp::self != bp::self );
        { //::SireVol::CoordGroupArray::operator=
        
            typedef ::SireVol::CoordGroupArray & ( ::SireVol::CoordGroupArray::*assign_function_type)( ::SireVol::CoordGroupArray const & ) ;
            assign_function_type assign_function_value( &::SireVol::CoordGroupArray::operator= );
            
            CoordGroupArray_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CoordGroupArray_exposer.def( bp::self == bp::self );
        { //::SireVol::CoordGroupArray::operator[]
        
            typedef ::SireVol::CoordGroup const & ( ::SireVol::CoordGroupArray::*__getitem___function_type)( ::quint32 ) const;
            __getitem___function_type __getitem___function_value( &::SireVol::CoordGroupArray::operator[] );
            
            CoordGroupArray_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireVol::CoordGroupArray::remove
        
            typedef void ( ::SireVol::CoordGroupArray::*remove_function_type)( ::quint32 ) ;
            remove_function_type remove_function_value( &::SireVol::CoordGroupArray::remove );
            
            CoordGroupArray_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("i") )
                , "Remove the ith CoordGroup from the array" );
        
        }
        { //::SireVol::CoordGroupArray::remove
        
            typedef void ( ::SireVol::CoordGroupArray::*remove_function_type)( ::quint32,int ) ;
            remove_function_type remove_function_value( &::SireVol::CoordGroupArray::remove );
            
            CoordGroupArray_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("i"), bp::arg("count") )
                , "Remove count CoordGroups from the array, starting with the ith CoordGroup" );
        
        }
        { //::SireVol::CoordGroupArray::remove
        
            typedef void ( ::SireVol::CoordGroupArray::*remove_function_type)( ::QVarLengthArray< unsigned int, 256 > const & ) ;
            remove_function_type remove_function_value( &::SireVol::CoordGroupArray::remove );
            
            CoordGroupArray_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("idxs") )
                , "Remove the specified CoordGroups from the array" );
        
        }
        { //::SireVol::CoordGroupArray::rotate
        
            typedef void ( ::SireVol::CoordGroupArray::*rotate_function_type)( ::SireMaths::Quaternion const &,::SireMaths::Vector const & ) ;
            rotate_function_type rotate_function_value( &::SireVol::CoordGroupArray::rotate );
            
            CoordGroupArray_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("quat"), bp::arg("point") )
                , "Rotate all of the coordinates in this array using the quaternion\nquat around the point point" );
        
        }
        { //::SireVol::CoordGroupArray::rotate
        
            typedef void ( ::SireVol::CoordGroupArray::*rotate_function_type)( ::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
            rotate_function_type rotate_function_value( &::SireVol::CoordGroupArray::rotate );
            
            CoordGroupArray_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("rotmat"), bp::arg("point") )
                , "Rotate all of coordinates in this array using the matrix rotmat\nabout the point point" );
        
        }
        { //::SireVol::CoordGroupArray::rotate
        
            typedef void ( ::SireVol::CoordGroupArray::*rotate_function_type)( ::quint32,::SireMaths::Quaternion const &,::SireMaths::Vector const & ) ;
            rotate_function_type rotate_function_value( &::SireVol::CoordGroupArray::rotate );
            
            CoordGroupArray_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("i"), bp::arg("quat"), bp::arg("point") )
                , "Rotate all of the coordinates in the CoordGroup at index i using\nthe quaternion quat about the point point\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireVol::CoordGroupArray::rotate
        
            typedef void ( ::SireVol::CoordGroupArray::*rotate_function_type)( ::quint32,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
            rotate_function_type rotate_function_value( &::SireVol::CoordGroupArray::rotate );
            
            CoordGroupArray_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("i"), bp::arg("rotmat"), bp::arg("point") )
                , "Rotate all of the coordinates in the CoordGroup at index i using\nthe matrix rotmat about the point point\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireVol::CoordGroupArray::size
        
            typedef int ( ::SireVol::CoordGroupArray::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireVol::CoordGroupArray::size );
            
            CoordGroupArray_exposer.def( 
                "size"
                , size_function_value
                , "Return the number of CoordGroups in this array" );
        
        }
        { //::SireVol::CoordGroupArray::toString
        
            typedef ::QString ( ::SireVol::CoordGroupArray::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireVol::CoordGroupArray::toString );
            
            CoordGroupArray_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this array" );
        
        }
        { //::SireVol::CoordGroupArray::transform
        
            typedef void ( ::SireVol::CoordGroupArray::*transform_function_type)( ::SireMaths::Transform const & ) ;
            transform_function_type transform_function_value( &::SireVol::CoordGroupArray::transform );
            
            CoordGroupArray_exposer.def( 
                "transform"
                , transform_function_value
                , ( bp::arg("t") )
                , "Transform all of coordinates in this array using the transformation t" );
        
        }
        { //::SireVol::CoordGroupArray::transform
        
            typedef void ( ::SireVol::CoordGroupArray::*transform_function_type)( ::quint32,::SireMaths::Transform const & ) ;
            transform_function_type transform_function_value( &::SireVol::CoordGroupArray::transform );
            
            CoordGroupArray_exposer.def( 
                "transform"
                , transform_function_value
                , ( bp::arg("i"), bp::arg("t") )
                , "Transform all of the coordinates in the CoordGroup at index i using\nthe transformation t\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireVol::CoordGroupArray::translate
        
            typedef void ( ::SireVol::CoordGroupArray::*translate_function_type)( ::SireMaths::Vector const & ) ;
            translate_function_type translate_function_value( &::SireVol::CoordGroupArray::translate );
            
            CoordGroupArray_exposer.def( 
                "translate"
                , translate_function_value
                , ( bp::arg("delta") )
                , "Translate all of the coordinates in this array by delta" );
        
        }
        { //::SireVol::CoordGroupArray::translate
        
            typedef void ( ::SireVol::CoordGroupArray::*translate_function_type)( ::quint32,::SireMaths::Vector const & ) ;
            translate_function_type translate_function_value( &::SireVol::CoordGroupArray::translate );
            
            CoordGroupArray_exposer.def( 
                "translate"
                , translate_function_value
                , ( bp::arg("i"), bp::arg("delta") )
                , "Translate all of the coordinates in the ith CoordGroup by delta" );
        
        }
        { //::SireVol::CoordGroupArray::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireVol::CoordGroupArray::typeName );
            
            CoordGroupArray_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireVol::CoordGroupArray::update
        
            typedef void ( ::SireVol::CoordGroupArray::*update_function_type)( ::quint32,::SireVol::CoordGroup const & ) ;
            update_function_type update_function_value( &::SireVol::CoordGroupArray::update );
            
            CoordGroupArray_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("i"), bp::arg("cgroup") )
                , "Update the CoordGroup at index i so that it is equal to cgroup. Note\nthat cgroup must contain the same number of coordinates as the existing\nCoordGroup at this index\nThrow: SireError::invalid_index\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireVol::CoordGroupArray::update
        
            typedef void ( ::SireVol::CoordGroupArray::*update_function_type)( ::quint32,::QVector< SireMaths::Vector > const & ) ;
            update_function_type update_function_value( &::SireVol::CoordGroupArray::update );
            
            CoordGroupArray_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("i"), bp::arg("coords") )
                , "Update the CoordGroup at index i so that it has coordinates coords\nThere must contain the same number of coordinates as the existing\nCoordGroup at this index\nThrow: SireError::invalid_index\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireVol::CoordGroupArray::update
        
            typedef void ( ::SireVol::CoordGroupArray::*update_function_type)( ::quint32,::SireMaths::Vector const *,int ) ;
            update_function_type update_function_value( &::SireVol::CoordGroupArray::update );
            
            CoordGroupArray_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("i"), bp::arg("coords"), bp::arg("ncoords") )
                , "Update the CoordGroup at index i so that it has coordinates coords\n(there are ncoords coordinates in this array)\nThere must contain the same number of coordinates as the existing\nCoordGroup at this index\nThrow: SireError::invalid_index\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireVol::CoordGroupArray::what
        
            typedef char const * ( ::SireVol::CoordGroupArray::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireVol::CoordGroupArray::what );
            
            CoordGroupArray_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        CoordGroupArray_exposer.staticmethod( "typeName" );
        CoordGroupArray_exposer.def( "__copy__", &__copy__);
        CoordGroupArray_exposer.def( "__deepcopy__", &__copy__);
        CoordGroupArray_exposer.def( "clone", &__copy__);
        CoordGroupArray_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireVol::CoordGroupArray >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CoordGroupArray_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireVol::CoordGroupArray >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CoordGroupArray_exposer.def_pickle(sire_pickle_suite< ::SireVol::CoordGroupArray >());
        CoordGroupArray_exposer.def( "__str__", &__str__< ::SireVol::CoordGroupArray > );
        CoordGroupArray_exposer.def( "__repr__", &__str__< ::SireVol::CoordGroupArray > );
        CoordGroupArray_exposer.def( "__len__", &__len_size< ::SireVol::CoordGroupArray > );
    }

}
