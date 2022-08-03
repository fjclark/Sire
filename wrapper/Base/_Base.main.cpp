// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License


#include "boost/python.hpp"

#include "boost/python/suite/indexing/vector_indexing_suite.hpp"

#include "Helpers/clone_const_reference.hpp"

#include "Array2DBase.pypp.hpp"

#include "Array2D_double_.pypp.hpp"

#include "ArrayProperty_QString_.pypp.hpp"

#include "ArrayProperty_double_.pypp.hpp"

#include "ArrayProperty_int_.pypp.hpp"

#include "BooleanProperty.pypp.hpp"

#include "CPUID.pypp.hpp"

#include "ChunkedVector_double_.pypp.hpp"

#include "CombineProperties.pypp.hpp"

#include "DoubleArrayProperty.pypp.hpp"

#include "FlopsMark.pypp.hpp"

#include "GeneralUnitArrayProperty.pypp.hpp"

#include "GeneralUnitProperty.pypp.hpp"

#include "Incremint.pypp.hpp"

#include "IntegerArrayProperty.pypp.hpp"

#include "LengthProperty.pypp.hpp"

#include "LinkToProperty.pypp.hpp"

#include "LowerCaseString.pypp.hpp"

#include "MajorMinorVersion.pypp.hpp"

#include "MemInfo.pypp.hpp"

#include "NoMangling.pypp.hpp"

#include "NullProperty.pypp.hpp"

#include "NumberProperty.pypp.hpp"

#include "PackedArray2D_DoubleArrayProperty.pypp.hpp"

#include "PackedArray2D_DoubleArrayProperty_Array.pypp.hpp"

#include "PackedArray2D_IntegerArrayProperty.pypp.hpp"

#include "PackedArray2D_IntegerArrayProperty_Array.pypp.hpp"

#include "PackedArray2D_PropertyList.pypp.hpp"

#include "PackedArray2D_PropertyList_Array.pypp.hpp"

#include "PackedArray2D_QString_.pypp.hpp"

#include "PackedArray2D_QString_Array.pypp.hpp"

#include "PackedArray2D_QVariant_.pypp.hpp"

#include "PackedArray2D_QVariant_Array.pypp.hpp"

#include "PackedArray2D_StringArrayProperty.pypp.hpp"

#include "PackedArray2D_StringArrayProperty_Array.pypp.hpp"

#include "PackedArray2D_double_.pypp.hpp"

#include "PackedArray2D_double_Array.pypp.hpp"

#include "PackedArray2D_int_.pypp.hpp"

#include "PackedArray2D_int_Array.pypp.hpp"

#include "Process.pypp.hpp"

#include "Properties.pypp.hpp"

#include "Property.pypp.hpp"

#include "PropertyList.pypp.hpp"

#include "PropertyMap.pypp.hpp"

#include "PropertyName.pypp.hpp"

#include "Range.pypp.hpp"

#include "SimpleRange.pypp.hpp"

#include "StringArrayProperty.pypp.hpp"

#include "StringMangler.pypp.hpp"

#include "StringProperty.pypp.hpp"

#include "TempDir.pypp.hpp"

#include "TimeProperty.pypp.hpp"

#include "TrigArray2DBase.pypp.hpp"

#include "TrigArray2D_double_.pypp.hpp"

#include "TrimString.pypp.hpp"

#include "UnitTest.pypp.hpp"

#include "UpperCaseString.pypp.hpp"

#include "VariantProperty.pypp.hpp"

#include "Version.pypp.hpp"

#include "_Base_free_functions.pypp.hpp"

#include "vector_less__double__greater_.pypp.hpp"

namespace bp = boost::python;

#include "SireBase_containers.h"

#include "SireBase_registrars.h"

#include "SireBase_properties.h"

#include "SireBase/propertymap.h"

#include "SireBase/stringproperty.h"

#include "SireBase/numberproperty.h"

#include "SireBase/lengthproperty.h"

#include "SireBase/propertylist.h"

#include <QString>

#include "SireBase/slice.h"

void autoconvert_Slice();

BOOST_PYTHON_MODULE(_Base){
    register_SireBase_objects();

    register_SireBase_containers();

    register_vector_less__double__greater__class();

    register_Array2DBase_class();

    register_Array2D_double__class();

    register_Property_class();

    register_ArrayProperty_QString__class();

    register_ArrayProperty_double__class();

    register_ArrayProperty_int__class();

    register_BooleanProperty_class();

    register_CPUID_class();

    register_ChunkedVector_double__class();

    register_CombineProperties_class();

    register_DoubleArrayProperty_class();

    register_FlopsMark_class();

    register_GeneralUnitArrayProperty_class();

    register_GeneralUnitProperty_class();

    register_Incremint_class();

    register_IntegerArrayProperty_class();

    register_LengthProperty_class();

    register_LinkToProperty_class();

    register_StringMangler_class();

    register_LowerCaseString_class();

    register_MajorMinorVersion_class();

    register_MemInfo_class();

    register_NoMangling_class();

    register_NullProperty_class();

    register_NumberProperty_class();

    register_PackedArray2D_QString__class();

    register_PackedArray2D_QVariant__class();

    register_PackedArray2D_DoubleArrayProperty_class();

    register_PackedArray2D_IntegerArrayProperty_class();

    register_PackedArray2D_PropertyList_class();

    register_PackedArray2D_StringArrayProperty_class();

    register_PackedArray2D_double__class();

    register_PackedArray2D_int__class();

    register_Process_class();

    register_Properties_class();

    register_PropertyList_class();

    register_PropertyMap_class();

    register_PropertyName_class();

    register_Range_class();

    register_SimpleRange_class();

    register_StringArrayProperty_class();

    register_StringProperty_class();

    register_TempDir_class();

    register_TimeProperty_class();

    register_TrigArray2DBase_class();

    register_TrigArray2D_double__class();

    register_TrimString_class();

    register_UnitTest_class();

    register_UpperCaseString_class();

    register_VariantProperty_class();

    register_Version_class();

    register_PackedArray2D_QString_Array_class();

    register_PackedArray2D_QVariant_Array_class();

    register_PackedArray2D_DoubleArrayProperty_Array_class();

    register_PackedArray2D_IntegerArrayProperty_Array_class();

    register_PackedArray2D_PropertyList_Array_class();

    register_PackedArray2D_StringArrayProperty_Array_class();

    register_PackedArray2D_double_Array_class();

    register_PackedArray2D_int_Array_class();

    register_SireBase_properties();

    bp::implicitly_convertible< QString, SireBase::PropertyName >();

    bp::implicitly_convertible< SireBase::Property, SireBase::PropertyName >();

    bp::implicitly_convertible< QHash<QString,SireBase::PropertyName>, SireBase::PropertyMap >();

    bp::implicitly_convertible< SireBase::DoubleArrayProperty, SireBase::StringArrayProperty >();

    bp::implicitly_convertible< SireBase::DoubleArrayProperty, SireBase::IntegerArrayProperty >();

    bp::implicitly_convertible< SireBase::DoubleArrayProperty, SireBase::PropertyList >();

    bp::implicitly_convertible< SireBase::IntegerArrayProperty, SireBase::StringArrayProperty >();

    bp::implicitly_convertible< SireBase::IntegerArrayProperty, SireBase::DoubleArrayProperty >();

    bp::implicitly_convertible< SireBase::IntegerArrayProperty, SireBase::PropertyList >();

    bp::implicitly_convertible< SireBase::StringArrayProperty, SireBase::PropertyList >();

    autoconvert_Slice();

    register_free_functions();
}

