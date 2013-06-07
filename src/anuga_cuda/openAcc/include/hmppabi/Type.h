/*
 * Copyright (C) 2008-2013 CAPS entreprise.  All Rights Reserved.
 * 
 * The source code contained or described herein and all documents  related
 * to the source code ("Material") are owned  by  CAPS  entreprise  or  its
 * suppliers or licensors.
 * 
 * Title to the Material remains with CAPS entreprise or its suppliers  and
 * licensors.  The Material contains trade secrets and proprietary and con-
 * fidential information of CAPS entreprise or its suppliers and licensors.
 * 
 * The Material is protected by the French intellectual property code,  in-
 * tellectual property laws and international treaties.  No part of the Ma-
 * terial may be used, copied, reproduced, modified,  published,  uploaded,
 * posted, transmitted, distributed or disclosed in any  way  without  CAPS
 * entreprise's prior express written permission.
 * 
 * No license under any patent, copyright, trade secret or other  intellec-
 * tual property right is granted to or conferred upon you by disclosure or
 * delivery of the Material, either expressly, by implication,  inducement,
 * estoppel or otherwise.
 * 
 * Any license under such intellectual property rights  must  be  expressed
 * and approved by CAPS entreprise in writing.
 */
#ifndef HMPPABI_TYPE_H
#define HMPPABI_TYPE_H

#include <string>
#include <vector>
#include <map>

#include <hmppabi/MemorySpace.h>
#include <hmppabi/Intent.h>

namespace hmppabi
{

enum Kind
{
KIND_NONE = 0,
#define D(num, cname) cname = num,
#include <hmppabi/KindType.def>
#undef D
KIND_VOID = 100,
KIND_MAX
};

enum Scalar
{
SCALAR_NONE = 0,
#define D(num, name, ctype, dtyte, h5type, ffi_type, signed_, unsigned_) name = num,
#include <hmppabi/ScalarType.def>
#undef D
SCALAR_MAX
};


enum DimensionKind { UNDEFINED_DIMENSION, INTEGER_DIMENSION, IDENTIFIER_DIMENSION  };

struct ShapeElement
{
  ShapeElement(); // construct an empty shape
  ShapeElement(const std::string & name); // construct a shape element which is a string
  ShapeElement(int digit); // construct a shape element which is an integer
  std::string unparse(); // unparse the shape element

  bool isInteger() { return kind == INTEGER_DIMENSION; }
  bool isIdentifier() { return kind == IDENTIFIER_DIMENSION; }

  DimensionKind kind;
  size_t integer;
  std::string identifier;
};

class Type
{
public:

public:
  Type();
  Type(Kind k);
  // FIXME: sam> est-ce que ce constructeur sous-entend que le champ scalar est utilisé avec différents kind possibles ?
  Type(Kind k, Scalar s);
  Type(Kind k, const Type & subtype);
  ~Type();

public:
  static Type from_dpil_text(const std::string & s);
  static Type from_c_text(const std::string & s);

  static Type make_pointer(int level, Type t0);
  static Type make_complex(Type t);
  static Type make_signed(Type t);
  static Type make_unsigned(Type t);

  std::string asDPILText();
  std::string asCText();

  Type flattenize();
  Type bind(const std::map<std::string, int> & symbols);

public:
  Kind getKind() const { return kind_; }
  void setKind(Kind kind);

  Scalar getScalar() const { return scalar_; }
  void setScalar(Scalar scalar);

  size_t getSize() const;

  Type & getPointedType();
  std::vector<Type> & getSubTypes() { return subtypes_; }
  void addSubType(const Type & type);

  std::vector<ShapeElement> & getShape() { return shape_; }
  unsigned int getNumElement() const;
  void setShape(const std::vector<ShapeElement> & shape) { shape_ = shape; }
  void addShapeElement(const ShapeElement & se);

  std::string getDeclareName() const { return declare_name_; }
  void setDeclareName(const std::string & name);
  bool hasDeclareName() const { return ! declare_name_.empty(); }

  std::string getTypeName() const;
  void setTypeName(const std::string & name);
  bool hasTypeName() const;

  void setMemorySpace(MemorySpace memspace);
  MemorySpace getMemorySpace() { return memory_space_; }

  void setIntent(Intent intent) { intent_ = intent; }
  Intent getIntent() { return intent_; }

  bool isComplex() const;

  bool isArray() const { return getKind()==KIND_ARRAY; }
  bool isPointer() const { return getKind()==KIND_POINTER; }
  bool isScalar() const { return getKind()==KIND_SCALAR; }
  bool isStruct() const { return getKind()==KIND_STRUCT; }

private:
  std::vector<Type> flattenize_type();

private:
  Kind kind_;
  Scalar scalar_;
  std::string declare_name_;
  std::string type_name_;
  MemorySpace memory_space_;
  Intent intent_;
  std::vector<Type> subtypes_;
  std::vector<ShapeElement> shape_;
};

} // namespace hmppabi


#endif






