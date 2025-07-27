/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016, 2017 Danny van Dyk
 *
 * Copied from the Paludis package manager, which is
 * Copyright (c) 2008-2013 Ciaran McCreesh
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef EOS_GUARD_EOS_UTILS_VISITOR_HH
#define EOS_GUARD_EOS_UTILS_VISITOR_HH 1

#include <eos/utils/visitor-fwd.hh>

#include <functional>

namespace eos
{
    template <typename Visitor_> class AcceptVisitor
    {
        private:
            Visitor_ & _v;

        public:
            using result_type = void;

            /// @name Visitor operations
            /// @{
            AcceptVisitor(Visitor_ & v) :
                _v(v)
            {
            }

            template <typename T_>
            void
            operator() (T_ & t) const
            {
                t.accept(_v);
            }

            /// @}
    };

    template <typename Visitor_, typename Returning_> class AcceptVisitorReturning
    {
        private:
            Visitor_ & _v;

        public:
            using result_type = Returning_;

            /// @name Visitor operations
            /// @{
            AcceptVisitorReturning(Visitor_ & v) :
                _v(v)
            {
            }

            template <typename T_>
            Returning_
            operator() (T_ & t) const
            {
                return t.template accept_returning<Returning_>(_v);
            }

            /// @}
    };

    /*!
     * Convenience function for using a visitor with a standard algorithm.
     */
    template <typename Visitor_>
    AcceptVisitor<Visitor_>
    accept_visitor(Visitor_ & v)
    {
        return AcceptVisitor<Visitor_>(v);
    }

    template <typename Returning_, typename Visitor_>
    AcceptVisitorReturning<Visitor_, Returning_>
    accept_visitor_returning(Visitor_ & v)
    {
        return AcceptVisitorReturning<Visitor_, Returning_>(v);
    }

    template <typename> struct ExtractFirstArgumentType;

    template <typename T_, typename R_, typename A1_, typename... As_> struct ExtractFirstArgumentType<R_ (T_::*)(A1_, As_...) const>
    {
            using Type = A1_;
    };

    template <typename T_> using FirstCallArgumentType = typename ExtractFirstArgumentType<decltype(&T_::operator())>::Type;

    template <typename> struct ExtractResultType;

    template <typename T_, typename R_, typename... As_> struct ExtractResultType<R_ (T_::*)(As_...) const>
    {
            using Type = R_;
    };

    template <typename T_> using CallResultType = typename std::remove_const<typename ExtractResultType<decltype(&T_::operator())>::Type>::type;

    template <typename Revisitor_, typename Result_, typename... Cases_> struct MadeVisitor;

    template <typename Revisitor_, typename Result_> struct MadeVisitor<Revisitor_, Result_>
    {
            Result_ visit(const NoType<0u> &) const;
    };

    template <typename> struct CallThisCaseNeedsTwoArgs;

    template <typename T_, typename R_, typename A1_> struct CallThisCaseNeedsTwoArgs<R_ (T_::*)(A1_) const>
    {
            enum
            {
                value = false
            };
    };

    template <typename T_, typename R_, typename A1_, typename A2_> struct CallThisCaseNeedsTwoArgs<R_ (T_::*)(A1_, A2_) const>
    {
            enum
            {
                value = true
            };
    };

    template <typename Result_, typename Case_, typename V_, bool needs_two_args_> struct CallThisCase;

    template <typename Result_, typename Case_, typename V_> struct CallThisCase<Result_, Case_, V_, false>
    {
            static Result_
            call(const Case_ & thiscase, const FirstCallArgumentType<Case_> & v, const V_ &)
            {
                return thiscase(v);
            }
    };

    template <typename Result_, typename Case_, typename V_> struct CallThisCase<Result_, Case_, V_, true>
    {
            static Result_
            call(const Case_ & thiscase, const FirstCallArgumentType<Case_> & v, const V_ & revisitor)
            {
                return thiscase(v, accept_visitor_returning<Result_>(revisitor));
            }
    };

    template <typename Case_, typename V_> struct CallThisCase<void, Case_, V_, true>
    {
            static void
            call(const Case_ & thiscase, const FirstCallArgumentType<Case_> & v, const V_ & revisitor)
            {
                thiscase(v, accept_visitor(revisitor));
            }
    };

    template <typename Revisitor_, typename Result_, typename Case_, typename... Rest_>
    struct MadeVisitor<Revisitor_, Result_, Case_, Rest_...> : MadeVisitor<Revisitor_, Result_, Rest_...>
    {
            const Case_ & thiscase;

            MadeVisitor(const Case_ & c, const Rest_ &... cases) :
                MadeVisitor<Revisitor_, Result_, Rest_...>(cases...),
                thiscase(c)
            {
            }

            using MadeVisitor<Revisitor_, Result_, Rest_...>::visit;

            Result_
            visit(const FirstCallArgumentType<Case_> & v) const
            {
                return CallThisCase<Result_, Case_, Revisitor_, CallThisCaseNeedsTwoArgs<decltype(&Case_::operator())>::value>::call(thiscase,
                                                                                                                                     v,
                                                                                                                                     *static_cast<const Revisitor_ *>(this));
            }
    };

    template <typename Result_, typename... Cases_> struct BaseMadeVisitor : MadeVisitor<BaseMadeVisitor<Result_, Cases_...>, Result_, Cases_...>
    {
            BaseMadeVisitor(const Cases_ &... cases) :
                MadeVisitor<BaseMadeVisitor<Result_, Cases_...>, Result_, Cases_...>(cases...)
            {
            }
    };

    template <typename Case_, typename... Cases_>
    auto
    make_visitor(const Case_ & firstcase, const Cases_ &... cases) -> BaseMadeVisitor<CallResultType<Case_>, Case_, Cases_...>
    {
        return BaseMadeVisitor<CallResultType<Case_>, Case_, Cases_...>{ firstcase, cases... };
    }

    template <typename Result_, typename Base_> using Revisit = std::function<Result_(const Base_ &)>;

    template <> class DeclareAbstractVisitMethods<TypeListTail>
    {
        public:
            void forward_visit(const NoType<0u> &);
    };

    template <typename TypeList_> class DeclareAbstractVisitMethods : public virtual DeclareAbstractVisitMethods<typename TypeList_::Tail>
    {
        public:
            using DeclareAbstractVisitMethods<typename TypeList_::Tail>::forward_visit;

            virtual void forward_visit(typename TypeList_::Item &) = 0;
    };

    template <typename TypeList_> class WrappedVisitorBase : public virtual DeclareAbstractVisitMethods<TypeList_>
    {};

    template <typename RealClass_> class ImplementVisitMethods<RealClass_, TypeListTail>
    {
        public:
            void forward_visit(const NoType<1u> &);
    };

    template <typename RealClass_, typename TypeList_>
    class ImplementVisitMethods : public virtual DeclareAbstractVisitMethods<TypeList_>, public ImplementVisitMethods<RealClass_, typename TypeList_::Tail>
    {
        public:
            using ImplementVisitMethods<RealClass_, typename TypeList_::Tail>::forward_visit;

            virtual void
            forward_visit(typename TypeList_::Item & n)
            {
                static_cast<RealClass_ *>(this)->perform_visit(n);
            }
    };

    template <typename TypeList_, typename UnwrappedVisitor_>
    class WrappedVoidResultVisitor : public WrappedVisitorBase<TypeList_>, public ImplementVisitMethods<WrappedVoidResultVisitor<TypeList_, UnwrappedVisitor_>, TypeList_>
    {
        private:
            UnwrappedVisitor_ & _unwrapped_visitor;

        public:
            WrappedVoidResultVisitor(UnwrappedVisitor_ & v) :
                _unwrapped_visitor(v)
            {
            }

            template <typename C_>
            void
            perform_visit(C_ & t)
            {
                _unwrapped_visitor.visit(t);
            }
    };

    template <typename TypeList_, typename Result_, typename UnwrappedVisitor_>
    class WrappedNonVoidResultVisitor :
        public WrappedVisitorBase<TypeList_>,
        public ImplementVisitMethods<WrappedNonVoidResultVisitor<TypeList_, Result_, UnwrappedVisitor_>, TypeList_>
    {
        private:
            UnwrappedVisitor_ & _unwrapped_visitor;

        public:
            Result_ result;

            WrappedNonVoidResultVisitor(UnwrappedVisitor_ & v, const Result_ & r) :
                _unwrapped_visitor(v),
                result(r)
            {
            }

            template <typename C_>
            void
            perform_visit(C_ & t)
            {
                result = _unwrapped_visitor.visit(t);
            }
    };

    template <typename BaseClass_, typename VisitableTypeList_> class DeclareAbstractAcceptMethods
    {
        private:
            virtual void _real_accept(WrappedVisitorBase<VisitableTypeList_> &)                                               = 0;
            virtual void _real_accept_const(WrappedVisitorBase<typename MakeTypeListConst<VisitableTypeList_>::Type> &) const = 0;

        public:
            using VisitableTypeList  = VisitableTypeList_;
            using VisitableBaseClass = BaseClass_;

            template <typename UnwrappedVisitor_>
            void
            accept(UnwrappedVisitor_ & v)
            {
                WrappedVoidResultVisitor<VisitableTypeList_, UnwrappedVisitor_> vv(v);
                _real_accept(vv);
            }

            template <typename UnwrappedVisitor_>
            void
            accept(UnwrappedVisitor_ & v) const
            {
                WrappedVoidResultVisitor<typename MakeTypeListConst<VisitableTypeList_>::Type, UnwrappedVisitor_> vv(v);
                _real_accept_const(vv);
            }

            template <typename UnwrappedVisitor_>
            void
            accept(const UnwrappedVisitor_ & v)
            {
                WrappedVoidResultVisitor<VisitableTypeList_, const UnwrappedVisitor_> vv(v);
                _real_accept(vv);
            }

            template <typename UnwrappedVisitor_>
            void
            accept(const UnwrappedVisitor_ & v) const
            {
                WrappedVoidResultVisitor<typename MakeTypeListConst<VisitableTypeList_>::Type, const UnwrappedVisitor_> vv(v);
                _real_accept_const(vv);
            }

            template <typename Result_, typename UnwrappedVisitor_>
            Result_
            accept_returning(UnwrappedVisitor_ & v, const Result_ & r = Result_())
            {
                WrappedNonVoidResultVisitor<VisitableTypeList_, Result_, UnwrappedVisitor_> vv(v, r);
                _real_accept(vv);
                return vv.result;
            }

            template <typename Result_, typename UnwrappedVisitor_>
            Result_
            accept_returning(const UnwrappedVisitor_ & v, const Result_ & r = Result_())
            {
                WrappedNonVoidResultVisitor<VisitableTypeList_, Result_, const UnwrappedVisitor_> vv(v, r);
                _real_accept(vv);
                return vv.result;
            }

            template <typename Result_, typename UnwrappedVisitor_>
            Result_
            accept_returning(UnwrappedVisitor_ & v, const Result_ & r = Result_()) const
            {
                WrappedNonVoidResultVisitor<typename MakeTypeListConst<VisitableTypeList_>::Type, Result_, UnwrappedVisitor_> vv(v, r);
                _real_accept_const(vv);
                return vv.result;
            }

            template <typename Result_, typename UnwrappedVisitor_>
            Result_
            accept_returning(const UnwrappedVisitor_ & v, const Result_ & r = Result_()) const
            {
                WrappedNonVoidResultVisitor<typename MakeTypeListConst<VisitableTypeList_>::Type, Result_, const UnwrappedVisitor_> vv(v, r);
                _real_accept_const(vv);
                return vv.result;
            }

            template <typename Case_, typename... Cases_>
            auto
            make_accept_returning(const Case_ & firstcase, const Cases_ &... cases) const -> CallResultType<Case_>
            {
                return this->accept_returning<CallResultType<Case_>>(make_visitor(firstcase, cases...));
            }

            template <typename... Cases_>
            void
            make_accept(const Cases_ &... cases) const
            {
                this->accept(make_visitor(cases...));
            }
    };

    template <typename BaseClass_, typename RealClass_>
    class ImplementAcceptMethods : public virtual DeclareAbstractAcceptMethods<BaseClass_, typename BaseClass_::VisitableTypeList>
    {
        private:
            void
            _real_accept(WrappedVisitorBase<typename BaseClass_::VisitableTypeList> & v)
            {
                v.forward_visit(*static_cast<RealClass_ *>(this));
            }

            void
            _real_accept_const(WrappedVisitorBase<typename MakeTypeListConst<typename BaseClass_::VisitableTypeList>::Type> & v) const
            {
                v.forward_visit(*static_cast<const RealClass_ *>(this));
            }
    };
} // namespace eos

#endif
