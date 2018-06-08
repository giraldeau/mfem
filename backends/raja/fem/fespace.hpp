// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability see http://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#ifndef MFEM_BACKENDS_RAJA_FE_SPACE_HPP
#define MFEM_BACKENDS_RAJA_FE_SPACE_HPP

#include "../../../config/config.hpp"
#if defined(MFEM_USE_BACKENDS) && defined(MFEM_USE_RAJA)

namespace mfem
{

namespace raja
{

/// TODO: doxygen
class RajaFiniteElementSpace : public mfem::PFiniteElementSpace
{
protected:

   Layout e_layout;
   int globalDofs, localDofs;
   int vdim;
   mfem::Ordering::Type ordering;
   raja::array<int> globalToLocalOffsets;
   raja::array<int> globalToLocalIndices,*reorderIndicess;
   raja::array<int> localToGlobalMap;
   mfem::Operator *restrictionOp, *prolongationOp;
public:
   /// TODO: doxygen
   RajaFiniteElementSpace(const Engine &e, mfem::FiniteElementSpace &fespace);

   /// Virtual destructor
   virtual ~RajaFiniteElementSpace();

   /// TODO: doxygen
   const Engine &RajaEngine() const
   { return *static_cast<const Engine *>(engine.Get()); }

   /// TODO: doxygen
   raja::device GetDevice(int idx = 0) const
   { return RajaEngine().GetDevice(idx); }

   mfem::Mesh* GetMesh() const { return fes->GetMesh(); }

   mfem::FiniteElementSpace* GetFESpace() const { return fes; }

   Layout &RajaVLayout() const
   { return *fes->GetVLayout().As<Layout>(); }

   Layout &RajaTrueVLayout() const
   { return *fes->GetTrueVLayout().As<Layout>(); }

   Layout &RajaEVLayout() { return e_layout; }

#ifdef MFEM_USE_MPI
   bool isDistributed() const { return (RajaEngine().GetComm() != MPI_COMM_NULL); }
#else
   bool isDistributed() const { return false; }
#endif

   bool hasTensorBasis() const
   { return dynamic_cast<const mfem::TensorBasisElement*>(fes->GetFE(0)); }

   mfem::Ordering::Type GetOrdering() const { return ordering; }

   int GetGlobalDofs() const { return globalDofs; }
   int GetLocalDofs() const { return localDofs; }

   int GetDim() const { return fes->GetMesh()->Dimension(); }
   int GetVDim() const { return vdim; }

   int GetVSize() const { return globalDofs * vdim; }
   int GetTrueVSize() const { return fes->GetTrueVSize(); }
   int GetGlobalVSize() const { return globalDofs*vdim; /* FIXME: MPI */ }
   int GetGlobalTrueVSize() const { return fes->GetTrueVSize(); }

   int GetNE() const { return fes->GetNE(); }

   const mfem::FiniteElementCollection* FEColl() const
   { return fes->FEColl(); }
   const mfem::FiniteElement* GetFE(const int idx) const
   { return fes->GetFE(idx); }

   //const int* GetElementDofMap() const { return elementDofMap; }
   //const int* GetElementDofMapInverse() const { return elementDofMapInverse; }

   const mfem::Operator* GetRestrictionOperator() { return restrictionOp; }
   const mfem::Operator* GetProlongationOperator() { return prolongationOp; }

   const raja::array<int> GetLocalToGlobalMap() const
   { return localToGlobalMap; }

   void GlobalToLocal(const raja::Vector &globalVec, Vector &localVec) const;
   void LocalToGlobal(const Vector &localVec, Vector &globalVec) const;

};
   
} // namespace mfem::raja

} // namespace mfem

#endif // defined(MFEM_USE_BACKENDS) && defined(MFEM_USE_RAJA)

#endif // MFEM_BACKENDS_RAJA_FE_SPACE_HPP
