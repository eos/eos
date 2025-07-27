/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_EOS_UTILS_VISITOR_FWD_HH
#define EOS_GUARD_EOS_UTILS_VISITOR_FWD_HH 1

#include <eos/utils/type-list-fwd.hh>

namespace eos
{
    template <typename TypeList_> class DeclareAbstractVisitMethods;

    template <> class DeclareAbstractVisitMethods<TypeListTail>;

    template <typename TypeList_> class WrappedVisitorBase;

    template <typename RealClass_, typename TypeList_> class ImplementVisitMethods;

    template <typename RealClass_> class ImplementVisitMethods<RealClass_, TypeListTail>;

    template <typename TypeList_, typename UnwrappedVisitor_> class WrappedVoidResultVisitor;

    template <typename TypeList_, typename Result_, typename UnwrappedVisitor_> class WrappedNonVoidResultVisitor;

    template <typename BaseClass_, typename VisitableTypeList_> class DeclareAbstractAcceptMethods;

    template <typename BaseClass_, typename RealClass_> class ImplementAcceptMethods;

    template <unsigned u_> class NoType;
} // namespace eos
#endif
